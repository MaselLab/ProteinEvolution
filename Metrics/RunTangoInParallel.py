import csv, os, sys, json, mysql.connector, time, datetime
from multiprocessing import Pool, Process
from Bio import SeqIO
from Bio.Seq import Seq

'''
Author : Sara Willis
Date   : February 11, 2019
--------------------------

The purpose of this script is to read in protein sequences from a MySQL database and to calculate aggregation metrics for them using the executable Tango. 

This script makes use of the multiprocessing module to run Tango on multiple proteins in the database simultaneously. This allows the script to be run on fairly large datatables in a reasonable period of time. 


--------------------------------------------------------------------------------------
Submodules

There are only two submodules that need to be imported from the submodule script (RunTangoInParallel_Submodules.py):

   1) LoadMySQLConnectionInformation -- Stores the data that is needed to access the MySQL database where the protein sequences are stored. There is no input for this function and there are four outputs:
           1) Username to access the MySQL database
           2) Password used to access the MySQL database
           3) IP address used to acces the MySQL database
           4) The name of the database that needs to be accessed

   2) CalculateAggregation_Tango -- This calculates the various aggregation metrics for the user and returns them to the user. The input/output is the following

      Input:
           1) ProteinSequence -- The input amino acid sequence 
           2) batch=0 -- If the process is being run in parallel, this is used to differentiate the various files that are output from Tango
           3) PfamIndices= None -- If the metrics are desired for domains, the start/stop indices of the domain should be included as a list ([start,stop]). Multiple domains are allowed.
           4) Domains = False -- This lets the program know whether only the full protein is being considered (False), or if domains are being considered (True). The default is False. If it is set to True, PfamIndices is required
           5) AggThreshold=5 -- This is the threshold for a raw Tango score to be considered aggregation prone
           6) RunThreshold=5 -- This is the number of scores above the AggThreshold needed in a row to be considered an APR
           7) ct='N' -- Required by Tango, c-terminal
           8) nt='N', Required by Tango, n-terminal
           9) ph='7.0' -- Required by Tango, temperature (K)
           10) te='298.15' -- Required by Tango, temperature (K)
           11) io='0.02' -- Reqired by Tango, ionic strength

     Output:
           1) AggregationList -- This is output as a list of the raw scores output by Tango. NULL if Tango was not run on this sequence. Either the sequence was too long (>2000 residues), there were more than three unknown residues in a row (X's) in the sequence, or an unknown residue was found in an APR.
           2) NumOfAggProneRegions_FullProtein -- The total number of APRs showing up in the sequence
           3) NumOfAggProneRegions_Density_FullProtein -- The total number of APRs divided by the length of the sequence
           4) NumOfAAInAPRs_FullProtein -- The number of amino acids inside an APR
           5) NumOfAAInAPRs_Density_FullProtein -- The number of amino acids inside an APR divided by the length of the sequence

--------------------------------------------------------------------------------------
           The output entries below are only returned in the input field Domains is set to True and PfamIndices are supplied. All output is supplied as a list with each entry matching the index in the input Pfam list

           6) numOfAggProneRegions_List_Pfam -- The number of APRs in the domain 
           7) numOfAggProneRegions_Density_List_Pfam -- The number of APRs in the domain divided by the length of the domain
           8) numOfAAsInAPRs_List_Pfam -- The number of residues in APRs in the domain
           9) numOfAAsInAPRs_Density_List_Pfam -- The number of residues in the APRs in the domain divided by the length of the domain
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
    numberOfProcesses = 5

    # User Paths
    # ------------------------
    

    # User's MySQL Database information
    # ------------------------
    Database = ''                                           # Name of user's database
    User = ''                                               # Username to access Fusion/MySQL
    Host = ''                                               # MySQL host
    Password = ''                                           # User's MySQL password
    DataTable = 'TestUploadTable_Protein_UseMeForTests'     # Name of the DataTable in the Database where sequences are stored
    ProteinColumnName = 'ProteinSequence'                   # Name of the column where the user's proteins are stored in the data table
    UIDColumnName = 'UID'                                   # Name of the column where the user's Table UIDs are stored in the data table (type = string or int)
    RawTangoScoresColumnName = 'RawTangoOutputScores'       # Name of the column where the raw tango scores should be saved (type = text)
    NumberOfAPRsColumnName = 'NumberOfAPRs'                 # Name of the column where the number of APRs should be saved (type = int)
    DensityOfAPRsColumnName = 'DensityOfAggProneRegions'    # Name of the column where the density of APRs should be saved (type = float)
    NumberOfAAsInAPRsColumnName = 'NumberOfAAsInAPRs'       # Name of the column where the number of amino acids in APRs should be saved (type = int)
    DensityOfAAsInAPRsColumnName = 'DensityOfAAsInAPRs'     # Name of the column where the density of amino acids in APRs should be saved (type = float)
    



    return {'NumberOfProcesses':numberOfProcesses,'Database':Database, 'User':User, 'Host':Host, 'Password':Password, 'DataTable':DataTable, 'ProteinColumnName': ProteinColumnName, 'UIDColumnName':UIDColumnName, 'RawTangoScoresColumnName':RawTangoScoresColumnName, 'NumberOfAPRsColumnName':NumberOfAPRsColumnName, 'DensityOfAPRsColumnName':DensityOfAPRsColumnName, 'NumberOfAAsInAPRsColumnName':NumberOfAAsInAPRsColumnName, 'DensityOfAAsInAPRsColumnName':DensityOfAAsInAPRsColumnName}



################################################################################
#                                   SUBMODULES                                 #
################################################################################




##### Partition MySQL Database #####

'''
This function partitions the datatable into roughly equal-sized chunks to be processed in parallel. The number of chunks is determined by the number of processes defined by the user in UserOptions
'''
def PartitionDatabase(DatabaseEntries):
    
    # The number of entries that need to be processed is determined
    numberOfEntries = len(DatabaseEntries)
    # The database will be partitioned into equal chunks saved as sublists in a master list
    partitionedResults = []
    # The size of each sublist is determined by divididing the number of entries in the full table by the number of processes
    sizeOfChunk =  numberOfEntries/UserOptions()['NumberOfProcesses']
    # The size of the chunk is then rounded to an integer so it can be used to index the list from the database
    sizeOfChunk = int(sizeOfChunk)
    # A list of indices to partition the list from the database is generated
    chunks = list(range(0,numberOfEntries,sizeOfChunk))
    # Since the final chunk may not be as long as the others, if the last index of the database list is not included in the partitioning, it is included
    if chunks[-1] < numberOfEntries:
        chunks.append(numberOfEntries)
        # Then, the portion of the database with indices lying between each chunk index is added to a list and saved in a master list
    for i in range(len(chunks)-1):
        start = chunks[i]
        stop = chunks[i+1]
        partitionedResults.append(list(results[start:stop]))
    # In order to run Tango and get the output, temporary files need to be created and read in. As a result, things go awry if the various batchs are not kept track of when Tango is run in parallel. As a result, a batch number is added to each sublist so that it can be used by Tango to differentiate filenames. See the submodules file to see where Batch comes in.
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
    # It will look for runs of X's (unknown residues) of three or more in length.
    XRun = 0    
    # Each residue in the sequence is searched
    for residue in Sequence:
        # If an unknown residue is found, a run begins
        if residue == 'X':
            # Each X in succession adds one to the length of the run
            XRun += 1            
        # If the residue is known, the "run" ends. This could either be a "true" run (where XRun =/= 0), or it could follow another known residue. Either way, each time a known residue is found, the run length is checked
        else:            
            # If the length of the run is greater than or equal to 3, the sequence does not pass this quality filter and the check is immediately returned as False for the user
            if XRun >= 3:
                return False           
            # So long as the run is below 3, the check continues and XRun is reset to 0
            else:
                XRun = 0        
    # If no runs >= 3 were found, the sequence passes the quality check and the user is returned True             
    return True


# ----------------------------------------------------------------------------------





##### Replace X's in Raw Tango Output #####

'''
The purpose of this submodule is to take in a file of raw Tango output and return to the user a list of the concatenated raw scores. In the process, it reintroduces unknown amino acids located in the body of the original sequence back into the raw score list as X's
'''

def ReplaceXTangoOutput(TangoFilename,ProteinSequence):
    
    # A header exists in the Tango output file, so this row is skipped
    Header = True    
    # We define an empty list where our raw Tango values will be stored
    AggregationList = []    
    # We use an index to keep track of where we are in the body of the protein sequence. We need to do this since the list of aggregation scores is less than or equal to the length of the protein sequence (since there's no score where an unknown amino acid exists), so we need to keep track of where we are.
    proteinIndex = 0    
    # We open the file
    with open(TangoFilename, 'r') as f:        
        # Define a reader to parse each row
        reader = csv.reader(f, delimiter = '\t')
        # Then search each line in the output file
        for row in reader:            
            # We skip the first heading row
            if Header == True:
                Header = False
            # Then begin to search through the raw data
            else:                
                # We define variable names for the raw output in each row
                ResidueNumber = row[0]
                AA = row[1]
                Beta = row[2]
                Turn= row[3]
                Helix = row[4]               
                # We're only interested in the Aggregation score for this analysis
                Aggregation = float(row[5])
                Conc_Stab_Aggregation = row[6]                
                # We check the aggregation score against the AA at our location in the body of the protein sequence. 
                if ProteinSequence[proteinIndex] != 'X':                    
                    # If it's not unknown (an X), then we add the aggregation score to our output list
                    AggregationList.append(Aggregation)                    
                    # And we add one to the protein index so that parsing the next row, we'll be at the appropriate location in the protein sequence
                    proteinIndex += 1                   
                # If the amino acid is unknown, then we need to make sure to add in an X to the list out output aggregation scores before we append the aggregation score. First, we check to make sure that there's only one X that precedes the aggregation score (there cannot be more than two since we excluded sequences with three or more sequential unknown sequences prior to our analyses).
                elif ProteinSequence[proteinIndex:proteinIndex+2] != 'XX':                    
                    # If only one X precedes our score, we add an X to our list of aggregation scores
                    AggregationList.append('X')                    
                    # We then add this row's aggregation score
                    AggregationList.append(Aggregation)
                    # and add 2 to the protein index so we're at the appropriate location when we look at the next row
                    proteinIndex += 2
                else: 
                    # If two XX's precede our aggregation score, we add both to the aggregation score list and then the current aggregation score
                    AggregationList += ['X','X',Aggregation]                 
                    # We then add 3 to our proteinIndex so we're at the appropriate place during our next iteration
                    proteinIndex += 3                    
        # Once we're done parsing the file, we make sure that we've included all the X's in our Aggregation (since if the protein sequence ends with X's, these will not be accounted for in our for loop)
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
    # Each aggregation score is checked in the raw list
    for aggValue in AggregationList:
        # a run is strung together from aggregation scores and X's
        if (type(aggValue) == float and float(aggValue) > AggThreshold) or aggValue == 'X':
            # The element in the run is added to the run list
            run.append(aggValue)
        # As soon as a score that doesn't exceed the AggTheshold is found, the run is broken
        else:
            # If the length of the run exceeds the RunThreshold, the run is checked to make sure there are no unknown AAs in the body of the run
            if len(run) >= RunThreshold:
                # Flanking unknown AAs are removed from the run and do not contribute to the length of the APR. As opposed to unknown AAs in the body of the run, these do not disqualify a sequence from analysis
                if run[0:2] == 'XX':
                    run = run[2:]
                if run[0] == 'X':
                    run = run[1:]
                if run[-2:] == 'XX':
                    run = run[:-2]
                if run[-1:] == 'X':
                    run = run[:-1]
                # After the flanking X's are removed from the run, if the length of the run exceeds the RunThreshold, then the sequence is checked to make sure there are no X's in the body of the sequence
                if len(run) >= RunThreshold:
                    if 'X' in run:
                        return False
                    else:
                        # If the run contains no unknown amino acids, then it is kept
                        runs.append(run)
                # If the run does not exceed the RunThreshold after the flanking unknown amino acids are removed, then the run is discarded and the aggregation list continues to be searched
                else:
                    pass
                runCount = 0
                run = []
    return runs


# ----------------------------------------------------------------------------------





##### Calculate Pfam Aggregation #####


def CalculatePfamAggregation(AggregationList, AggregationRuns, PfamIndices):
    # The number of aggregation prone regions, number of AAs in APRs, and the density of both will be found for each Pfam region
    numOfAAsInAPRs_List = []
    numOfAAsInAPRs_Density_List = []
    numOfAggProneRegions_List = []
    numOfAggProneRegions_Density_List = []
    
    for pfam in PfamIndices:
        # The number of AAs located in APRs is initially set to zero for each pfam region
        numOfAAsInAPRs = 0
        # Similarly, the number of agg-prone regions is set to zero at the beginning of each pfam run
        numberOfRegions = 0
        # Each run of Agg-prone AAs is then looked at to see whether it overlaps with the pfam region
        for run in AggregationRuns:
            lengthOfRun = len(run)
            # The APR is located within the body of the full aggregation list
            for ind in (i for i,e in enumerate(AggregationList) if e == run[0]):
                if AggregationList[ind:ind+lengthOfRun]:
                    startIndex = ind
                    stopIndex = ind + lengthOfRun
            # The number of AAs located within the APR that intersect with the Pfam region is found
            Intersection = set(range(pfam[0],pfam[1]+1)).intersection(set(range(startIndex, stopIndex+1)))
            # If the APR doesn't overlap with the Pfam domain at all, it is ignored
            if len(Intersection) == 0:
                pass
            # If the APR intersects with the Pfam domain, the number of AAs in that APR that are also located within the Pfam region are noted, and the number of Aggregation-prone regions in that Pfam region goes up by one
            else:
                numOfAAsInAPRs += len(Intersection)
                numberOfRegions += 1
        # Once each one of the APRs has been searched, the density values are found
        numOfAAsInAPRs_Density = numOfAAsInAPRs/(pfam[1]-pfam[0])
        numOfAggProneRegions_Density =numberOfRegions/(pfam[1]-pfam[0])
    
        # Then, because we need scores for all of the Pfam regions, we add the scores to a list. They will be added in the order that the pfam IDs appear in so they can be coordinated with the correct region
        numOfAAsInAPRs_List.append(numOfAAsInAPRs)
        numOfAAsInAPRs_Density_List.append(numOfAAsInAPRs_Density)
        numOfAggProneRegions_List.append(numberOfRegions)
        numOfAggProneRegions_Density_List.append(numOfAggProneRegions_Density)
    # Once each pfam region is searched, the various lists are returned to the user
    return numOfAAsInAPRs_List, numOfAAsInAPRs_Density_List, numOfAggProneRegions_List, numOfAggProneRegions_Density_List


# ----------------------------------------------------------------------------------


##### Calculate Pfam Aggregation #####

'''
This is a simple submodule that returns null values to the user in the event that Tango fails for some reason. This function can be called and immediately returned to the user. This is done to reduce clutter in the main script (writing this all would only lead to a mess in the main body, so it's defined here).

It takes two entries for input 

   1) Domain -- this is either True or False, and is declared in the main function that it's called by. This determines the number of entries returned to the user (there will be no list of Null values returned corresponding to Pfam regions if Domain = False)

   2) PfamIndices -- This is either a list of tuples corresponding to the locations of Pfam IDs within the body of the protein, or None. Either way, it's defined by the user in the main program 
'''

def ReturnNullValues(Domain, PfamIndices):
    if Domain == True:
        AggregationList = 'NULL'
        NumAggProneRegions_FullProtein = 'NULL'
        NumOfAggProneRegions_Density_FullProtein = 'NULL'
        NumOfAAInAPRs_FullProtein = 'NULL'
        NumOfAAInAPRs_Density_FullProtein = 'NULL'
        numOfAggProneRegions_List_Pfam = ['NULL' for i in PfamIndices]
        NumOfAggProneRegions_Density_List_Pfam = ['NULL' for i in PfamIndices]
        numOfAAsInAPRs_List_Pfam = ['NULL' for i in PfamIndices]
        numOfAAsInAPRs_Density_List_Pfam = ['NULL' for i in PfamIndices]
        return AggregationList, NumAggProneRegions_FullProtein, NumOfAggProneRegions_Density_FullProtein, NumOfAAInAPRs_FullProtein, NumOfAAInAPRs_Density_FullProtein, numOfAggProneRegions_List_Pfam, NumOfAggProneRegions_Density_List_Pfam, numOfAAsInAPRs_List_Pfam, numOfAAsInAPRs_Density_List_Pfam

    else:
        AggregationList = 'NULL'
        NumAggProneRegions_FullProtein = 'NULL'
        NumOfAggProneRegions_Density_FullProtein = 'NULL'
        NumOfAAInAPRs_FullProtein = 'NULL'
        NumOfAAInAPRs_Density_FullProtein = 'NULL'
        return AggregationList, NumAggProneRegions_FullProtein, NumOfAggProneRegions_Density_FullProtein, NumOfAAInAPRs_FullProtein, NumOfAAInAPRs_Density_FullProtein


# ----------------------------------------------------------------------------------



##### Upload Tango Aggregation Metrics To MySQL #####
'''

This function runs the 
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
    # n will keep track of the number of proteins processed. 
    n = 0
    # The batch number is extracted 
    batch = results[-1]
    # and then the batch number is removed from the list
    results = results[:-1]
    #  each UID, protein pair in the list is then read in
    for result in results:
        # The number of proteins processed goes up by one
        n+=1
        # Each 1000 entries processed notifies the user of the number processed so far and how long it took to process that number
        if (n%10000) == 0:
            print('Number Of Entries Processed: %s\nTime Taken: %s     '%(n,datetime.datetime.now()-start_time),end='\r')
        # The UID is then extracted from the tuple
        UID = result[1]
        # As well as the Protein sequence
        ProteinSequence = result[0]
        
        # And Tango is run on the Protein sequence. For a description of the inputs/outputs of this function, see the top of this script
        AggregationList, NumOfAggProneRegions, NumOfAggProneRegions_Density, NumOfAAsInAPRs,NumOfAAsInAPRs_Density = CalculateAggregation_Tango(ProteinSequence,batch)
        # Sometimes the Aggregation list is returned as the string NULL for a variety of reasons (see top of this script for these reasons). If this is not the case, a comma-separated string is defined consisting of all the raw scores output from Tango. If it's NULL, then it's left alone.
        if AggregationList != 'NULL':
            AggregationList = ','.join([str(i) for i in AggregationList])
        # And the MySQL table is updated
        # Because MySQL is finicky about dealing with syntax in update/insert statements, two insert statements are defined, one assuming the UID in the user's
        # data table is an integer. If this fails, the syntax is changed so the insertion statement assumes the user's UID is a string
        try:
            # First the script assumes the UID is an integer
            UpdateStatement = "UPDATE %s SET %s = '%s', %s='%s', %s='%s',%s='%s',%s='%s' WHERE %s=%s" %(UserOptions()['DataTable'],UserOptions()['RawTangoScoresColumnName'],AggregationList,UserOptions()['NumberOfAPRsColumnName'],NumOfAggProneRegions,UserOptions()['DensityOfAPRsColumnName'],NumOfAggProneRegions_Density,UserOptions()['NumberOfAAsInAPRsColumnName'],NumOfAAsInAPRs,UserOptions()['DensityOfAAsInAPRsColumnName'],NumOfAAsInAPRs_Density,UserOptions()['UIDColumnName'],UID)
            mycursor.execute(UpdateStatement)
        except:
            # If the insertion statement fails, the script assumes the UID is a string
            UpdateStatement = "UPDATE %s SET %s = '%s', %s='%s', %s='%s',%s='%s',%s='%s' WHERE %s='%s'" %(UserOptions()['DataTable'],UserOptions()['RawTangoScoresColumnName'],AggregationList,UserOptions()['NumberOfAPRsColumnName'],NumOfAggProneRegions,UserOptions()['DensityOfAPRsColumnName'],NumOfAggProneRegions_Density,UserOptions()['NumberOfAAsInAPRsColumnName'],NumOfAAsInAPRs,UserOptions()['DensityOfAAsInAPRsColumnName'],NumOfAAsInAPRs_Density,UserOptions()['UIDColumnName'],UID)
            mycursor.execute(UpdateStatement)
        # The results are then committed to the mysql table
        cnx.commit()
    # Once all the proteins have been processed, the connection is closed
    cnx.close()
    

# ----------------------------------------------------------------------------------




'''
The purpose of this script is to calculate the aggregation propensity of a protein sequence, as well as of defined Pfam polypeptide subsequences. 


'''

# The function is defined. It takes two mandatory parameters as input and includes seven optional parameters that are defined with default parameters if not defined by the user
def CalculateAggregation_Tango(ProteinSequence, batch=0,PfamIndices= None, Domains = False,AggThreshold=5, RunThreshold=5,ct='N',nt='N',ph='7.0', te='298.15', io='0.02'):
    
    # The output from the Tango run is redirected to a flat file so that it can be parsed 
    OutputFilename = 'TangoOutput_%s.txt'%batch

    # The protein sequences are transformed to get rid of unknown amino acids so that they can be processed by Tango
    ProteinSequence = ProteinSequence.replace('U','X').replace('B','X').replace('Z','X').replace('O','X').replace('J','X')
    ProteinSequence_NoUnknownResidues = ProteinSequence.replace('X','').replace(' ','')

    # If a sequence has three or more unknown amino acids in a row, then it is excluded from processing by Tango. Similarly, sequences with more than 2000 amino acids are also excluded because they cannot be processed by Tango. The sequence that has to pass this second filter is the sequence once all unknown residues are removed. This is because this is the actual sequence that gets fed into Tango
    if SequencePassed_TangoQualityControl(ProteinSequence) == False or len(ProteinSequence_NoUnknownResidues) >= 2000:

        # Null values will be returned to the user using a function which alters the null values that are returned based on whether only a full protein is being run through the function, or a protein with pfam regions (as defined by the input variable Domain)
        return ReturnNullValues(Domains, PfamIndices)
        
    # Otherwise, if the sequence passes the first round of filters, it is fed into Tango using the command-line interface. The data are saved to the OutputFilename. The command-line output is saved to a temporary file so that it can be deleted each round. This is to prevent too much being printed to the command line for large queries
    os.system('./tango_x86_64_release %s ct=N nt=N ph=7.0 te=298.15 io=0.02 seq=%s >redirectedOutput_%s.txt' %(OutputFilename.replace('.txt',''),ProteinSequence_NoUnknownResidues,batch))

    # If, for whatever reason, output from Tango is not generated, then null values are returned to the user and the function exits
    if os.path.exists(OutputFilename) == False:
        print('Error running Tango, no output file found.\nReturning null values')
        return ReturnNullValues(Domains, PfamIndices)

    # If the redirected stdout file exists, it's removed
    if os.path.exists('redirectedOutput_%s.txt'%batch):
        os.remove('redirectedOutput_%s.txt'%batch)

    # The Tango raw data are then processed and the X's are returned to the raw list where they appear in the protein sequence. The user is returned the list of raw Tango scores
    AggregationList = ReplaceXTangoOutput(OutputFilename,ProteinSequence)
    
    # The aggregation prone regions are then found and a list of lists is returned to the user, each sublist containing the values of each amino acid in the aggregation-prone run
    AggregationProneRegions = FindAPRs(AggregationList, AggThreshold, RunThreshold)
    if AggregationProneRegions == False:
        return ReturnNullValues(Domains, PfamIndices)

    # The number of aggregation prone regions in the protein is then the number of runs in the AggregationProneRegions list
    NumOfAggProneRegions_FullProtein = len(AggregationProneRegions)
    
    # And the density is calculated as the number of aggregation prone regions divided by the length of the protein
    NumOfAggProneRegions_Density_FullProtein = NumOfAggProneRegions_FullProtein/len(AggregationList)
    
    # The number of amino acids in aggregatio-prone regions is then found using the list of runs. The length of each run is found and the lengths are summed
    NumOfAAInAPRs_FullProtein = sum(map(len,AggregationProneRegions))
    
    # And the density of AAs in APRs is found by dividing the total number of AAs in APRs by the length of the protein
    NumOfAAInAPRs_Density_FullProtein = NumOfAAInAPRs_FullProtein/len(AggregationList)

    # Because running Pfam regions through the pipeline is optional, only if Domains is set to True in the input will the program find the number of AAs in APRs, the number of APRs and the various density values for Pfam regions
    if Domains == True:
        
        # The values are calculated for the Pfam regions using a function and returned to the user
        numOfAAsInAPRs_List_Pfam, numOfAAsInAPRs_Density_List_Pfam, numOfAggProneRegions_List_Pfam, numOfAggProneRegions_Density_List_Pfam = CalculatePfamAggregation(AggregationList, AggregationProneRegions, PfamIndices)
        
        # and all values are returned to the user
        return AggregationList, NumOfAggProneRegions_FullProtein, NumOfAggProneRegions_Density_FullProtein, NumOfAAInAPRs_FullProtein, NumOfAAInAPRs_Density_FullProtein, numOfAggProneRegions_List_Pfam, numOfAggProneRegions_Density_List_Pfam, numOfAAsInAPRs_List_Pfam, numOfAAsInAPRs_Density_List_Pfam
    
    else:
        # If pfam regions were not looked for, only the data for the full regions is returned to the user
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
