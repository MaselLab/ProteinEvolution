import csv, os, sys, json, mysql.connector, time, datetime, random
from multiprocessing import Pool, Process
from Bio import SeqIO
from Bio.Seq import Seq

'''
Author : Sara Willis
Date   : March 11, 2019
--------------------------

The purpose of this script is to read in protein sequences from a MySQL database, to generate randomly scrambled peptide sequences from those proteins, and to calculate aggregation metrics for the randomized sequences using the executable Tango. This script was created with the intention of generating values to determine delta tango for protein sequences.

The scrambled peptides will be saved with their Tango output scores in a MySQL data table along with the UID that corresponds to the protein the randomized peptide was generated from. The user has the option to generate more than one randomized peptide sequence per protein using the NumberOfScrambledSequences variable in UserOptions. The table where the scrambled peptides and their metrics may or may not exist prior to running this script. If the user wants the script to generate the scrambled peptide table for them, set CreateRandomizedPeptideTable to True.

This script makes use of the multiprocessing module to run Tango on multiple proteins in the database simultaneously. This allows the script to be run on fairly large datatables in a reasonable period of time. 

Dependencies
------------

The user will need the executable Tango in order for this script to run: http://tango.crg.es/

The user will also need two non-standard Python modules that can be installed with the conda command:
   - BioPython       : https://anaconda.org/anaconda/biopython
   - mysql.connector : https://anaconda.org/anaconda/mysql-connector-python

'''


####################################################################################
#                                  User Options                                    #
####################################################################################

'''
This should be the only part of the script that the user needs to modify

The user options are stored as a function so that they can be used with ease throughout the rest of the script
'''

def UserOptions():
    # Number of Processes - Used for running script in parallel on large databases. Suggested maximum: 20
    # ------------------------
    numberOfProcesses = 20

    # User Paths
    # ------------------------
    TangoPath = './tango_x86_64_release'
    
    # User's MySQL Database information
    # ------------------------
    Database = ''                                           # Name of user's database
    User = ''                                               # Username to access Fusion/MySQL
    Host = ''                                               # MySQL host
    Password = ''                                           # User's MySQL password
    
    DataTable = 'Genomes_Multicellular_Protein'             # Name of the DataTable where protein sequences are stored
    ProteinColumnName = 'ProteinSequence'                   # Name of the column where the user's proteins are stored
    UIDColumnName = 'UID'                                   # Name of the UID column
    
    NumberOfScrambledSequences = 1                          # Define the number of randomly scrambled sequences the user wishes to generate per protein

    ScrambledTable = 'Genomes_Multicellular_ScrambledSequencesWithTango'  # Name of the column where the scrambled protein sequences will be stored with the Tango data
    ScrambledProteinColumnName = 'ScrambledSequence'        # Name of the column where the user's scrambled peptides are stored
    ScrambledProteinUIDColumnName = 'ProteinTableUID'       # There should be an indexed column in the scrambled table where the protein table's UIDs will be stored
    RawTangoScoresColumnName = 'RawTangoOutputScores'       # Raw tango scores column (type = text)
    NumberOfAPRsColumnName = 'NumOfAPRs'                    # Number of APRs column (type = int)
    DensityOfAPRsColumnName = 'DensityOfAggProneRegions'    # Density of APRs column (type = decimal)
    NumberOfAAsInAPRsColumnName = 'NumberOfAAsInAPRs'       # Number of amino acids in APRs (type = int)
    DensityOfAAsInAPRsColumnName = 'DensityOfAAsInAPRs'     # Density of amino acids in APRs (type = decimal)

    CreateRandomizedPeptideTable = True                     # If this is set to True, the MySQL ScrambledTable will be generated for the user. If a table with the
                                                            # designated name already exists, the program will exit with a warning 

    Verbose = True                                          # How chatty the program should be
    



    return {'NumberOfProcesses' : numberOfProcesses,
            
            'TangoPath' : TangoPath,
            
            'Database' : Database,
            'User' : User,
            'Host' : Host,
            'Password' : Password,
            
            'DataTable' : DataTable,
            'ProteinColumnName' : ProteinColumnName,
            'UIDColumnName' : UIDColumnName,

            'NumberOfScrambledSequences' : NumberOfScrambledSequences,

            'ScrambledTable' : ScrambledTable,
            'ScrambledProteinColumnName' : ScrambledProteinColumnName,
            'ScrambledProteinUIDColumnName' : ScrambledProteinUIDColumnName,
            'RawTangoScoresColumnName' : RawTangoScoresColumnName,
            'NumberOfAPRsColumnName' : NumberOfAPRsColumnName,
            'DensityOfAPRsColumnName' : DensityOfAPRsColumnName,
            'NumberOfAAsInAPRsColumnName' : NumberOfAAsInAPRsColumnName,
            'DensityOfAAsInAPRsColumnName' : DensityOfAAsInAPRsColumnName,

            'CreateRandomizedPeptideTable' : CreateRandomizedPeptideTable,

            'Verbose' : Verbose}



################################################################################
#                                   SUBMODULES                                 #
################################################################################

##### Create Randomized Peptide Table #####
'''
The user has the option to let this script create the randomized peptide table for them if it doesn't already exist. The table will be create with the name specified as ScrambledTable in UserOptions and will be generated with the following format:

ScrambledTable Columns:

UID                           -- Primary key, auto-incremented large integer. Not permitted to be NULL
ScrambledProteinUIDColumnName -- Indexed column containing the UIDs from the source protein table, big int
ScrambledProteinColumnName    -- Where the scrambled peptide sequences will be stored, text
RawTangoScoresColumnName      -- text
NumberOfAPRs                  -- Integer, max permitted length = 5
DensityOfAPRsColumnName       -- Decimal, length = 5, decimals = 4
NumberOfAAsInAPRsColumnName   -- Integer, max permitted length = 5
DensityOfAPRsColumnName       -- Decimal, length = 5, demicals = 4

If a table with the specified name already exists, the program will quit with an error

'''
def CreateRandomizedPeptideTable():
    cnx = mysql.connector.connect(user = UserOptions()['User'],
                                  password = UserOptions()['Password'],
                                  host = UserOptions()['Host'],
                                  database = UserOptions()['Database'])
    # And defines a cursor so the tables can be interacted with
    mycursor = cnx.cursor(buffered = True)
    try:
        # To make sure no existing tables get overwritten, we check to make sure the table doesn't already exist
        mycursor.execute("SELECT 1 FROM "+UserOptions()['ScrambledTable']+" LIMIT 1")
        print('\n\nWARNING!\n\nYou are trying to create a table that already exists!\nCheck your table names and try again or set CreateRandomizedPeptideTable to False\n\n')
    except:
        CreateTableStatement = "CREATE TABLE "+UserOptions()['ScrambledTable'] +"(UID BIGINT(100) not null auto_increment primary key, "+UserOptions()['ScrambledProteinUIDColumnName']+" BIGINT(100), "+UserOptions()['ScrambledProteinColumnName']+" TEXT, "+UserOptions()['RawTangoScoresColumnName']+" TEXT, "+UserOptions()['NumberOfAPRsColumnName']+" INT(5), " +UserOptions()['DensityOfAPRsColumnName']+" DECIMAL(5,4), " + UserOptions()['NumberOfAAsInAPRsColumnName']+ " INT(5), " + UserOptions()['DensityOfAAsInAPRsColumnName']+ " DECIMAL(5,4), INDEX protein_table_UID ("+UserOptions()['ScrambledProteinUIDColumnName']+"))"
        mycursor.execute(CreateTableStatement)
        cnx.commit()
        cnx.close()
        return
    sys.exit(0)
    

# ----------------------------------------------------------------------------------




##### Partition MySQL Database #####

'''
This function partitions the datatable into roughly equal-sized chunks to be processed in parallel. The number of chunks is determined by the number of processes defined by the user in UserOptions
'''
def PartitionDatabase(DatabaseEntries):
    
    numberOfEntries = len(DatabaseEntries)
    # The database will be partitioned into equal chunks saved as sublists in a master list
    partitionedResults = []
    # The size of each sublist is determined by divididing the number of entries in the full table by the number of processes
    sizeOfChunk =  numberOfEntries/UserOptions()['NumberOfProcesses']
    # The size of the chunk is then rounded to an integer so it can be used to index the list from the database
    sizeOfChunk = int(sizeOfChunk)
    # A list of indices to partition the list from the database is generated
    chunks = list(range(0,numberOfEntries,sizeOfChunk))
    # It's possible the last entry won't be included in the final chunk, so we make sure we include it
    if chunks[-1] < numberOfEntries:
        chunks[-1] = numberOfEntries
    # Then, the portion of the database with indices lying between each chunk index is added to a list and saved in a master list
    for i in range(len(chunks)-1):
        start = chunks[i]
        stop = chunks[i+1]
        partitionedResults.append(list(results[start:stop]))
    # In order to run Tango and get the output, temporary files need to be created and read in. As a result, things go awry if the
    # various batchs are not kept track of when Tango is run in parallel. To circumvent this problem, a batch number is added to each
    # sublist so that it can be used by Tango to differentiate filenames
    batch = 0
    for element in partitionedResults:
        element.append(batch)
        batch += 1
    return partitionedResults


# ----------------------------------------------------------------------------------



##### Quality Control For Tango - Discard Sequences with X Triplicates #####

'''
This function is used to check to see whether a sequence passes imposed quality filters before it's run through Tango. Specifically, if a protein has three or more unknown residues in sequence, then the protein does not pass the filter and the user is notified to not use the sequence for Tango analysis. 

There is only one input: The protein sequence being considered for analysis

And there is only one output: True or False depending on whether the sequence passes the quality check or not.
'''

# The function takes in a protein sequence as input
def SequencePassed_TangoQualityControl(Sequence):
    if Sequence.count('XXX') >= 1:
        return False
    else:
        return True

# ----------------------------------------------------------------------------------



##### Replace X's in Raw Tango Output #####

'''
The purpose of this submodule is to take in a file of raw Tango output and return to the user a list of the concatenated raw scores. In the process, it reintroduces unknown amino acids located in the body of the original sequence back into the raw score list as X's
'''

def ReplaceXTangoOutput(TangoFilename,ProteinSequence):
    
    Header = True    
    # We define an empty list where our raw Tango values will be stored
    AggregationList = []    
    # We use an index to keep track of where we are in the body of the protein sequence. We need to do this since the list of
    # aggregation scores is less than or equal to the length of the protein sequence since there's no score where an unknown
    # amino acid exists.
    proteinIndex = 0
    
    with open(TangoFilename, 'r') as f:        
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:            
            if Header == True:
                Header = False
            else:                
                # We're only interested in the Aggregation score for this analysis
                Aggregation = float(row[5])
                # We check the aggregation score against the AA at our location in the body of the protein sequence. 
                if ProteinSequence[proteinIndex] != 'X':                    
                    # If it's not unknown, then we add the aggregation score to our output list
                    AggregationList.append(Aggregation)                    
                    # And we add one to the protein index so that parsing the next row, we'll be at the appropriate location
                    # in the protein sequence
                    proteinIndex += 1                   
                # If the amino acid is unknown, then we need to make sure to add in an X to the list out output aggregation scores
                # before we append the aggregation score. First, we check to make sure that there's only one X that precedes the
                # aggregation score (there cannot be more than two since we excluded sequences with three or more sequential
                # unknown sequences prior to our analyses).
                elif ProteinSequence[proteinIndex:proteinIndex+2] != 'XX':                    
                    # If only one X precedes our score, we add an X to our list of aggregation scores
                    AggregationList.append('X')                    
                    # We then add this row's aggregation score
                    AggregationList.append(Aggregation)
                    # and add 2 to the protein index so we're at the appropriate location when we look at the next row
                    proteinIndex += 2
                else: 
                    # If two XX's precede our aggregation score, we add both to the aggregation score list and then the current
                    # aggregation score
                    AggregationList += ['X','X',Aggregation]                 
                    # We then add 3 to our proteinIndex so we're at the appropriate place during our next iteration
                    proteinIndex += 3                    
        # Once we're done parsing the file, we make sure that we've included all the X's in our Aggregation (since if the protein
        # sequence ends with X's, these will not be accounted for in our for loop)
        # If the length of the aggregation list matches the length of the protein sequence, then nothing needs to be done
        if len(ProteinSequence) - len(AggregationList) == 0:
            pass        
        # If the aggregation list is one less than the protein sequence, then it ended with an X. As a result, an X is added to the aggregation list
        elif len(ProteinSequence) - len(AggregationList) == 1:
            AggregationList.append('X')           
        # If two X's exist at the end of the protein sequence, they are added to the aggregation list
        elif len(ProteinSequence) - len(AggregationList) == 2:
            AggregationList += ['X','X']           
    # Once the file has been parsed, it is removed
    os.remove(TangoFilename)   
    # And the list of raw aggregation scores with X's appended is returned to the user
    return AggregationList


# ----------------------------------------------------------------------------------


##### Find APRs #####

'''
The purpose of this submodule is to search through a list of raw Tango scores to find runs of aggregation-prone amino acids. The function calculates the density of the AAs in APRs in the sequence to the length of the sequence, as well as the number of AAs in APRs in the sequence. These two values are returned to the user as well as a list of the specific runs. This is in case the user needs to calculate the aggregation scores for subsequences of the protein as well.

***** This differs from Scott's script:

If an unknown amino acid is found within an APR, the sequence is excluded from analyses and the user is returned Null values. If unknown amino acids exist anywhere else in the 
'''

def FindAPRs(AggregationList, AggThreshold, RunThreshold): 
    run = []
    runs = []
    for aggValue in AggregationList:
        # a run is strung together from aggregation scores and X's
        if (type(aggValue) == float and float(aggValue) > AggThreshold) or aggValue == 'X':
            run.append(aggValue)
        # As soon as a score that doesn't exceed the AggTheshold is found, the run is broken
        else:
            # If the length of the run exceeds the RunThreshold, the run is checked to make
            # sure there are no unknown AAs in the body of the run
            if len(run) >= RunThreshold:
                # Flanking unknown AAs are removed from the run and do not contribute to the
                # length of the APR. As opposed to unknown AAs in the body of the run, these
                # do not disqualify a sequence from analysis
                if run[0:2] == 'XX':
                    run = run[2:]
                if run[0] == 'X':
                    run = run[1:]
                if run[-2:] == 'XX':
                    run = run[:-2]
                if run[-1:] == 'X':
                    run = run[:-1]
                # After the flanking X's are removed from the run, if the length of the run exceeds
                # the RunThreshold, then the sequence is checked to make sure there are no X's in
                # the body of the sequence
                if len(run) >= RunThreshold:
                    if 'X' in run:
                        return False
                    else:
                        # If the run contains no unknown amino acids, then it is kept
                        runs.append(run)
                # If the run does not exceed the RunThreshold after the flanking unknown amino acids
                # are removed, then the run is discarded and the aggregation list continues to be searched
                else:
                    pass
                runCount = 0
                run = []
    return runs




# ----------------------------------------------------------------------------------


##### Calculate Pfam Aggregation #####

'''
This is a simple submodule that returns null values to the user in the event that Tango fails for some reason. This function can be called and immediately returned to the user. This is done to reduce clutter in the main script (writing this all would only lead to a mess in the main body, so it's defined here).

'''

def ReturnNullValues():

    AggregationList = None
    NumAggProneRegions_FullProtein = None
    NumOfAggProneRegions_Density_FullProtein = None
    NumOfAAInAPRs_FullProtein = None
    NumOfAAInAPRs_Density_FullProtein = None
    return AggregationList, NumAggProneRegions_FullProtein, NumOfAggProneRegions_Density_FullProtein, NumOfAAInAPRs_FullProtein, NumOfAAInAPRs_Density_FullProtein


# ----------------------------------------------------------------------------------



##### Upload Tango Aggregation Metrics To MySQL #####
'''
This is the main function that's run in parallel and which makes use of all the other submodules. 
'''

# A function is defined that takes in a list of UIDs and the proteins they point to as tuples, plus a batch number as input.
# The function then runs Tango on each of the input proteins and uploads the aggregation metrics to the user's MySQL table by using the supplied UID to update the row
def UploadAggregationMetricsToMYSQL_Tango(results):
    # The function connects to the mysql database
    cnx = mysql.connector.connect(user = UserOptions()['User'],
                                  password = UserOptions()['Password'],
                                  host = UserOptions()['Host'],
                                  database = UserOptions()['Database'])
    # And defines a cursor so the tables can be interacted with
    mycursor = cnx.cursor(buffered = True)
    
    sequencesProcessed = 0
    # The batch number is extracted 
    batch = results[-1]
    # and then the batch number is removed from the list
    results = results[:-1]
    #  each UID, protein pair in the list is then read in
    for result in results:
        sequencesProcessed+=1
        # Each 1000 entries processed notifies the user of the number processed so far and how long it took to process that number
        # (if Verbose is set to True)
        if (sequencesProcessed%10000) == 0 and UserOptions()['Verbose'] == True:
            print('Number Of Entries Processed: %s\nTime Taken: %s     '%(sequencesProcessed,datetime.datetime.now()-start_time),end='\r')
            sys.stdout.flush()

        UID = result[1]
        ProteinSequence = result[0]
        # We'll only generate scrambled sequences from proteins that passed a quality check. This is to reduce unnecessary processing since
        # if the initial protein didn't pass the quality check, it does not have a Tango value associated with it so there's nothing to
        # compare the scrambled version to
        if ProteinSequence != None and SequencePassed_TangoQualityControl(ProteinSequence) == True and len(ProteinSequence) < 2000:
            # We'll generate as many randomly scrambled sequences as are defined by the user in UserOptions
            for i in range(0,UserOptions()['NumberOfScrambledSequences']):
                RandomSequenceGenerated = False
                # When we generate a scrambled sequence, we make sure that it passes the quality filters required to run Tango on it
                # if it doesn't, we generate a new one and continue this process until we get a scrambled sequence of sufficient quality
                while RandomSequenceGenerated == False:
                    ScrambledProteinSequence = ''.join(random.sample(ProteinSequence,len(ProteinSequence)))
                    RandomSequenceGenerated = SequencePassed_TangoQualityControl(ScrambledProteinSequence)

                # And Tango is run on the Protein sequence.
                AggregationList, NumOfAggProneRegions, NumOfAggProneRegions_Density, NumOfAAsInAPRs,NumOfAAsInAPRs_Density = CalculateAggregation_Tango(ScrambledProteinSequence,batch)
                # Sometimes the Aggregation list is returned as None for a variety of reasons. If this is not the case,
                # a comma-separated string is defined consisting of all the raw scores output from Tango. If it's = None then it's left alone.
                if AggregationList != None:
                    AggregationList = ','.join([str(i) for i in AggregationList])

                # We then insert the values into the scrambled table
                InsertionStatement = "INSERT INTO "+UserOptions()['ScrambledTable']+" ("+','.join([UserOptions()['ScrambledProteinUIDColumnName'],UserOptions()['ScrambledProteinColumnName'],UserOptions()['RawTangoScoresColumnName'],UserOptions()['NumberOfAPRsColumnName'],UserOptions()['DensityOfAPRsColumnName'],UserOptions()['NumberOfAAsInAPRsColumnName'],UserOptions()['DensityOfAAsInAPRsColumnName']])+") VALUES (%s,%s,%s,%s,%s,%s,%s)"
                InsertionData = (UID,ScrambledProteinSequence,AggregationList,NumOfAggProneRegions,NumOfAggProneRegions_Density,NumOfAAsInAPRs,NumOfAAsInAPRs_Density)
                mycursor.execute(InsertionStatement,InsertionData)

                cnx.commit()
    # Once all the proteins have been processed, the connection is closed
    cnx.close()


# ----------------------------------------------------------------------------------




'''
This function runs Tango and makes use of other, smaller functions to calculate the tango metrics associated with a peptide sequence
'''

def CalculateAggregation_Tango(ProteinSequence, batch=0, AggThreshold=5, RunThreshold=5, ct='N', nt='N', ph='7.0', te='298.15', io='0.02'):
    
    # The output from the Tango run is redirected to a flat file so that it can be parsed 
    OutputFilename = 'TangoOutput_%s.txt'%batch

    # The protein sequences are transformed to get rid of unknown amino acids so that they can be processed by Tango
    ProteinSequence = ProteinSequence.replace('U','X').replace('B','X').replace('Z','X').replace('O','X').replace('J','X')
    ProteinSequence_NoUnknownResidues = ProteinSequence.replace('X','').replace(' ','')

        
    # The sequence is fed into Tango using the command-line interface. The data are saved to the OutputFilename.
    # The command-line output is saved to a temporary file so that it can be deleted each round.
    #This is to prevent too much being printed to the command line for large queries
    os.system('%s %s ct=N nt=N ph=7.0 te=298.15 io=0.02 seq=%s >redirectedOutput_%s.txt' %(UserOptions()['TangoPath'],OutputFilename.replace('.txt',''),ProteinSequence_NoUnknownResidues,batch))

    # If, for whatever reason, output from Tango is not generated, then null values are returned to the user and the
    # function exits
    if os.path.exists(OutputFilename) == False:
        if UserOptions()['Verbose'] == True:
            print('Error running Tango, no output file found.\nReturning null values')
            os.stdout.flush()
        return ReturnNullValues()

    # If the redirected stdout file exists, it's removed
    if os.path.exists('redirectedOutput_%s.txt'%batch):
        os.remove('redirectedOutput_%s.txt'%batch)

    # The Tango raw data are then processed and the X's are returned to the raw list where they appear in the protein sequence.
    # The user is returned the list of raw Tango scores
    AggregationList = ReplaceXTangoOutput(OutputFilename,ProteinSequence)
    
    # The aggregation prone regions are then found and a list of lists is returned to the user, each sublist containing the values of each amino acid in the aggregation-prone run
    AggregationProneRegions = FindAPRs(AggregationList, AggThreshold, RunThreshold)
    if AggregationProneRegions == False:
        return ReturnNullValues()

    # The number of aggregation prone regions in the protein is then the number of runs in the AggregationProneRegions list
    NumOfAggProneRegions_FullProtein = len(AggregationProneRegions)
    
    # And the density is calculated as the number of aggregation prone regions divided by the length of the protein
    NumOfAggProneRegions_Density_FullProtein = NumOfAggProneRegions_FullProtein/len(AggregationList)
    
    # The number of amino acids in aggregation-prone regions is then found using the list of runs.
    # The length of each run is found and the lengths are summed
    NumOfAAInAPRs_FullProtein = sum(map(len,AggregationProneRegions))
    
    # And the density of AAs in APRs is found by dividing the total number of AAs in APRs by the length of the protein
    NumOfAAInAPRs_Density_FullProtein = NumOfAAInAPRs_FullProtein/len(AggregationList)


    return AggregationList, NumOfAggProneRegions_FullProtein, NumOfAggProneRegions_Density_FullProtein, NumOfAAInAPRs_FullProtein, NumOfAAInAPRs_Density_FullProtein

################################################################################









####################################################################################
############################# Program Executes Below ###############################
####################################################################################


'''
The actual execution of this script takes place below. It utilizes all the functions defined in the Submodules section as well as the user options. 

'''

# The time the script starts running
start_time = datetime.datetime.now()
print('Beginning Aggregation Analysis (Tango)')
print('Current Time: %s\n' %start_time)

# If the User wants to create the the scrambled peptide table, the program makes the attempt. The
# function CreateRandomizedPeptideTable will exit if the table name already exits
if UserOptions()['CreateRandomizedPeptideTable'] == True:
    if UserOptions()['Verbose'] == True:
        print('\nCreating Scrambled Peptide Table\n')
        sys.stdout.flush()
    CreateRandomizedPeptideTable()

# The script connects to the MySQL database using the connection information supplied by the user
cnx = mysql.connector.connect(user = UserOptions()['User'],
                              password = UserOptions()['Password'],
                              host = UserOptions()['Host'],
                              database = UserOptions()['Database'])
# And a cursor is defined so the user may interact with the datatables
mycursor = cnx.cursor(buffered = True)

# All protein sequences and their UIDs are pulled from the datatable
PullProteinsStatement = 'SELECT %s,%s FROM %s'%(UserOptions()['ProteinColumnName'],UserOptions()['UIDColumnName'],UserOptions()['DataTable'])
mycursor.execute(PullProteinsStatement)

# The results are saved as a list
results = mycursor.fetchall()
# And the connection to the database is closed
cnx.close()

# The database will be partitioned into equal chunks saved as sublists in a master list
partitionedResults = PartitionDatabase(results)


# Here we run the UploadAggregationMetricsToMYSQL_Tango function in parallel with the pool number being the user-supplied value numberOfProcesses
if __name__ == '__main__':
    p = Pool(UserOptions()['NumberOfProcesses'])
    p.map(UploadAggregationMetricsToMYSQL_Tango,partitionedResults)

# Once the script has completed running, the user is notified and the total time taken to run is displayed
print('\n\nScript Complete')
print('Time Taken: %s'%(datetime.datetime.now()-start_time))
