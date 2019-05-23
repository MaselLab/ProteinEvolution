import os, sys, json, csv, mysql.connector, datetime
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool

'''
This script is used to run IUPred2 on proteins stored in a MySQL database 

Dependencies:
--------------

--> anaconda3 (specifically python3) : https://www.anaconda.com/download/

--> mysql.connector : https://pynative.com/install-mysql-connector-python/

--> Biopython : https://anaconda.org/conda-forge/biopython


'''

####################################################################################
################################## User Options ####################################
####################################################################################

'''
This should be the only piece of the code that needs to be modified by the user. 

The user options are stored as a function so that they can be used with ease throughout the rest of the script
'''

def UserOptions():
    # Number of Processes - Used for running script in parallel on large databases. Suggested maximum: 20
    # ------------------------
    numberOfProcesses = 20

    # User Paths
    # ------------------------
    pathToPython3 = '~/anaconda3/bin/python3'               # This is the path to the python3 executable. IUPred2 requires Python3 to run
    IupredScript = './iupred2a/iupred2a.py'                 # The path from this script's directory to the IUPred2 python script

    # User Options for IUPred
    # ------------------------
    DisorderType = 'long'                                   # Disorder options for IUPred. Can either be long, short, or glob
    CysteinesIncluded = True                                # True = include cysteines in input protein for analysis; False = excise cysteines
    RunAnchor = False                                       # Iupred2 has the new option to run Anchor on protein sequences. True run Anchor and returns the raw output to the user

    # User's MySQL Database information
    # ------------------------
    Database = ''                                           # Name of user's database
    User = ''                                               # Username to access Fusion/MySQL
    Host = ''                                               # MySQL host
    Password = ''                                           # User's MySQL password
    DataTable = 'NCBIGenomes_ScrambledSequences_Complete'   # Name of the DataTable in the Database where sequences are stored
    ProteinColumnName = 'ScrambledSequence'                 # Name of the column where the user's proteins are stored in the data table
    UIDColumnName = 'UID'                                   # Name of the column where the user's table uids are stored in the data table
    IupredAverageColumnName = 'MeanISD_IUPred2_WithCys'     # Name of the column where the user's IUPred averages are stored
    IupredRawScoreColumnName = 'ISD_RawIUPredScores_CysIncluded'   # Name of the column where the user's raw IUPred scores are stored
    AnchorRawScoreColumnName = ''                           # Name of the column where the user's raw Anchor scores are stored. If anchor is not run, this can be left blank

    # The user-defined options are returned as a dictionary for use throughout the rest of the script
    return {'NumberOfProcesses': numberOfProcesses, 'PathToPython3': pathToPython3, 'IUPredScript': IupredScript, 'DisorderType': DisorderType, 'CysteinesIncluded' : CysteinesIncluded, 'RunAnchor' : RunAnchor, 'Database' : Database, 'User' : User, 'Host' : Host, 'Password' : Password, 'DataTable' : DataTable, 'ProteinColumnName' : ProteinColumnName, 'UIDColumnName' : UIDColumnName, 'IUPredAverageColumnName' : IupredAverageColumnName, 'IUPredRawScoreColumnName' : IupredRawScoreColumnName, 'AnchorRawScoreColumnName' : AnchorRawScoreColumnName}



####################################################################################
#################################### Submodules ####################################
####################################################################################


##### Create Fasta File #####

'''
This submodule is used to format a sequence as a fasta file. This is necessary for IUPred to be able to read in protein sequences to analyze

Input:
--------

   1) UID -- The unique ID associated with the protein sequence
   2) Protein -- The protein sequence 
   3) FastaFilename -- The filename the user wants the fasta file stored as
   4) (Optional) Description - Description of the protein sequence

Output:
--------

   1) A fasta file is saved for the user using the FastaFilename for the name/location of the file
'''

def CreateFastaFile(UID, Protein,FastaFilename, Description='UnknownDescription'):
    
    # we save sys.stdout as a new variable so that we can redirect what's printed
    # to the terminal to a file, and then revert to printing to the terminal at the
    # end of the function's use
    orig_stdout = sys.stdout
    # BioPython is used to generate a sequence record so that it's module SeqIO
    # can be called to write in fasta format
    record = SeqRecord(Seq(Protein), id= '%s'%UID, description = Description)
    # stdout is redirected to a file with the name chosen by the user
    sys.stdout = open('%s'%FastaFilename, 'a')
    print(record.format('fasta'))
    # Once it has been printed, sys.stdout is closed and the program will start
    # printing to the terminal again
    sys.stdout.close()
    # Future output is then redirected to the terminal
    sys.stdout=orig_stdout

    
# ----------------------------------------------------------------------------------

    
##### Calculate Disorder - Iupred 2 #####

'''
This is the submodule that performs the majority of the work to execute IUPred and to process the output

Input:
--------

   1) The UID of the sequence being processed
   2) The protein sequence to analyze

Output:
--------
   1) A comma-delimited string of all the raw IUPred disorder scores
   2) The average value of all raw IUPred disorder scores
   3) A comma-delimited string of all the raw Anchor scores (returns None if Anchor was not directed to run)
'''

def CalculateDisorder_IUPred2(UID,ProteinSequence, batch):

    # The cysteines are removed from the peptide sequence prior to IUPred analysis if the user sets the optional argument CysteinesIncluded to False
    if UserOptions()['CysteinesIncluded'] == False:
        ProteinSequence = ProteinSequence.replace('C','')
    # All stop codon characters are removed from the peptide sequence prior to analysis
    ProteinSequence = ProteinSequence.replace('*','')
    # And a temporary fasta filename is defined
    TempFastaFile = 'Temp_%s.fa'%batch
    # The filename for the output from IUPred is defined
    IupredOutputFile = 'IUPredOutputFile_%s.txt'%batch
    # Before any analysis takes place, it's first checked to make sure a temporary fasta file doesn't already exist. If it does, it's removed
    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)
    # A fasta file is created containing the input sequence with it's UID using the temporary fasta filename we defined
    CreateFastaFile(UID,ProteinSequence,TempFastaFile)
    # If any previous output files from an IUPred run exist, they are removed
    if os.path.exists(IupredOutputFile) == True:
        os.remove(IupredOutputFile)
    if UserOptions()['RunAnchor'] == False:
        # If the user doesn't want to run Anchor, then Iupred is executed normally using the options
        # specified at the top of this script
        os.system('%s %s %s %s >%s' %(UserOptions()['PathToPython3'],UserOptions()['IUPredScript'],TempFastaFile,UserOptions()['DisorderType'],IupredOutputFile))
    else:
        # If the user does want to run Anchor, then the user option "-a" is used during the execution of IUPred2
        os.system('%s %s -a %s %s >%s' %(UserOptions()['PathToPython3'],UserOptions()['IUPredScript'],TempFastaFile,UserOptions()['DisorderType'],IupredOutputFile))
    # Lists will be kept to keep track of the ISD/Anchor raw scores
    ISD_RawScore = []
    Anchor_RawScore = []
    # Once Iupred is run, the output file is opened and parsed to return the relevant information to the user
    with open(IupredOutputFile,'r') as f:
        # Be default, the output from IUPred is a tab-delimited file
        reader = csv.reader(f, delimiter = '\t')
        # Each row in the output file is searched
        for row in reader:
            # If the row starts with the character #, it is a human-readable row containing no data and so is skipped
            if len(row) in [0,1] or '#' in row[0]:
                pass
            else:
                # If the row doesn't start with #, then it contains data and so the row is parsed
                AANumber = int(row[0])                   # The first value in the row is the number of the amino acid in the peptide
                AA = row[1]                              # The peptide analyzed is the second value in the row
                ISDValue = row[2]                        # The ISD value for the particular AA is the third value in the row
                ISD_RawScore.append(str(format(float(ISDValue),'.4f')))            # We add the ISD value to the list of raw scores to be returned to the user
                # If the user has specified that Anchor should be run, a fourth value is pulled from the output file associated with Anchor
                if UserOptions()['RunAnchor'] == True:             
                    AnchorValue = row[3]                 # The Anchor value is the optional fourth element in the row.
                    Anchor_RawScore.append(AnchorValue)  # The anchor value is added to the raw anchor output list
    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)
    if os.path.exists(IupredOutputFile):
        os.remove(IupredOutputFile)
    # If Anchor is set to false, then only the disorder predictions of IUPred are returned to the user. The output is specifically:
    #    1) A comma-separated string of the raw output scores from IUPred
    #    2) The mean ISD calculated by averaging all raw output IUPred scores
    if UserOptions()['RunAnchor'] == False:
        return(','.join(ISD_RawScore),np.mean([float(i) for i in ISD_RawScore]),None)
    # If Anchor was specified to run, then a third element is returned to the user -- specifically a comma-separated string of the raw anchor output values
    else:
        return(','.join(ISD_RawScore),np.mean([float(i) for i in ISD_RawScore]),','.join(Anchor_RawScore))
    

# ----------------------------------------------------------------------------------
    

##### Partition Database #####

'''
The purpose of this submodule is to partition the database that IUPred is being run on into equal-sized chunks. The number of chunks is determined 
by the number of processes that the user has specified they want to utilize while running this script

Input: 
--------

   1) A list containing the data from each row in the user's database as a tuple

Output:
--------

   2) A master list containing the database partitioned into equal-sized chunks stored as sublists in the master list
'''

# First, a list of every entry in the database is fed in as input
def PartitionDataset(ProteinSequences):
    # And the total size of the database is determined
    numberOfSequences = len(ProteinSequences)
    # The partitioned database will be returned to the user as a list of lists
    partitionedSequences = []
    # The size of the chunk is an integer that's rounded from the total length of the database divided by the number
    # of processes the user wishes to use in the program's execution. Each process will be run on one of the sublists
    sizeOfChunk = numberOfSequences/UserOptions()['NumberOfProcesses']
    sizeOfChunk = int(sizeOfChunk)
    # The indices of the chunks are then determined by creating a list from 0 to the length of the database in steps
    # the size of the chunk
    chunks = list(range(0, numberOfSequences, sizeOfChunk))
    # If the final index is missing, then it's added to the chunk indices list
    if chunks[-1] < numberOfSequences:
        chunks.append(numberOfSequences)
    # Then, the database is partitioned using these indices
    for i in range(len(chunks)-1):
        start = chunks[i]
        stop = chunks[i+1]
        partitionedSequences.append(list(ProteinSequences[start:stop]))
    # Finally, a batch number is added to each of the sublists so that multiprocessing can be used
    batchNumber = 0
    for sublist in partitionedSequences:
        sublist.append(batchNumber)
        batchNumber+=1 
    # And the list of the partitioned database is returned to the user
    return partitionedSequences


# ----------------------------------------------------------------------------------


##### Run IUPred and Upload Results #####

'''
This submodule runs IUPred2 on a list of protein sequences using the user-specified options and uploads the results to the user's database.

Input:
------

    A list of protein sequences and their UIDs stored as tuples and with a final entry that specifies the batch number so this function can be run in parallel. 
    The batch number, in particular, is used to differentiate between output filenames. 

    The input, specifically, should have the following form:
    [ (UID_1,Protein_1), ... , (UID_N,Protein_N), batch_number ]

Output:
------

    Nothing is returned to the user. The results from this analysis are uploaded directly to the user's MySQL database using the options specified at the beginning of this script
'''
def RunIUPredAndUploadResults(ProteinSequences):
    start_time = datetime.datetime.now()

    # The batch number is removed from the input list and stored
    batch = ProteinSequences[-1]
    ProteinSequences =ProteinSequences[:-1]
    numberOfSequences = len(ProteinSequences)
    # A connection is then established with the MySQL server
    cnx = mysql.connector.connect(user = UserOptions()['User'],
                                  password = UserOptions()['Password'],
                                  host = UserOptions()['Host'],
                                  database = UserOptions()['Database'])
    # And a cursor is defined so the data tables can be interacted with 
    mycursor = cnx.cursor(buffered = True)
    # Each entry in the input list is then parsed and analyzed
    for ProteinSequence in ProteinSequences:
        # The first element in the tuple is the UID of the protein
        UID = ProteinSequence[0]
        # And the protein sequence is the second element
        protein = ProteinSequence[1]
        # And IUPred is run on the sequence
        ISDRawScores, IUPredMeanISD, AnchorRawScores = CalculateDisorder_IUPred2(UID, protein, batch)
        if UserOptions()['RunAnchor'] == True:
            # If the user has designated that they wish to run Anchor, the Anchor results are uploaded to the user's data table
            try:
                # mysql.connector is finicky about whether things are integers or strings.
                # First, the script assumes the user is using an integer for their UID in their datatable. If that raises an
                # exception, then the script attempts the upload statement assuming the UID is a string 
                UpdateSQLTableStatement = "UPDATE %s SET %s='%s',%s='%s',%s='%s' WHERE %s=%s"%(UserOptions()['DataTable'],UserOptions()['IUPredAverageColumnName'],IUPredMeanISD,UserOptions()['IUPredRawScoreColumnName'],ISDRawScores,UserOptions()['AnchorRawScoreColumnName'],AnchorRawScores,UserOptions()['UIDColumnName'],UID)
                mycursor.execute(UpdateSQLTableStatement)

            except:
                UpdateSQLTableStatement = "UPDATE %s SET %s='%s',%s='%s',%s='%s' WHERE %s='%s'"%(UserOptions()['DataTable'],UserOptions()['IUPredAverageColumnName'],IUPredMeanISD,UserOptions()['IUPredRawScoreColumnName'],ISDRawScores,UserOptions()['AnchorRawScoreColumnName'],AnchorRawScores,UserOptions()['UIDColumnName'],UID)
                mycursor.execute(UpdateSQLTableStatement)

        else:
            # Otherwise, only disorder predictions will be uploaded
            try:
                UpdateSQLTableStatement = "UPDATE %s SET %s='%s',%s='%s' WHERE %s=%s"%(UserOptions()['DataTable'],UserOptions()['IUPredAverageColumnName'],IUPredMeanISD,UserOptions()['IUPredRawScoreColumnName'],ISDRawScores,UserOptions()['UIDColumnName'],UID)
                mycursor.execute(UpdateSQLTableStatement)

            except:
                UpdateSQLTableStatement = "UPDATE %s SET %s='%s',%s='%s' WHERE %s='%s'"%(UserOptions()['DataTable'],UserOptions()['IUPredAverageColumnName'],IUPredMeanISD,UserOptions()['IUPredRawScoreColumnName'],ISDRawScores,UserOptions()['UIDColumnName'],UID)
                mycursor.execute(UpdateSQLTableStatement)

        
        # And commit the changes to the data table
        cnx.commit()
    # Once every protein has been processed, we close the connection
    cnx.close()
    print('Analysis Complete\nNumber of sequences processed: %s\nTime taken: %s'%(numberOfSequences,(datetime.datetime.now()-start_time)))
    
# ----------------------------------------------------------------------------------



####################################################################################
############################# Program Executes Below ###############################
####################################################################################

start_time = datetime.datetime.now()

print('Beginning IUPred2 Analysis\nNumber of Processes Being Used: %s\nCurrent Time: %s'%(UserOptions()['NumberOfProcesses'],start_time))

# We start by establishing a connection with the MySQL server

cnx = mysql.connector.connect(user = UserOptions()['User'],
                              password = UserOptions()['Password'],
                              host = UserOptions()['Host'],
                              database = UserOptions()['Database'])

# We define a cursor so we can interact with the data tables
mycursor = cnx.cursor(buffered = True)

# and we define an extraction statement to pull the relevant data from 
PullProteinsStatement = "SELECT %s,%s FROM %s" %(UserOptions()['UIDColumnName'],UserOptions()['ProteinColumnName'],UserOptions()['DataTable'])
# We then execute the command to extract the sequences
mycursor.execute(PullProteinsStatement)
# We pull all the results from the cursor using the fetchall command
ProteinSequences = mycursor.fetchall()
# And close the connection to the SQL database
cnx.close()
# And partition the database into equal sized chunks based on the number of processes the user
# wants to use
ProteinSequences = PartitionDataset(ProteinSequences)

# Finally, we use multiprocessing to run the script in parallel as many
# times as the user specifies
if __name__ == '__main__':
    # In the event that the partitioning of the database results in a list that has fewer sublists
    # than the number of processes the user wants, we modify the pool number to match the length
    # of partitioned database to avoid any errors
    if len(ProteinSequences) == UserOptions()['NumberOfProcesses']:
        p = Pool(UserOptions()['NumberOfProcesses'])
        p.map(RunIUPredAndUploadResults,ProteinSequences)
    else:
        p = Pool(len(ProteinSequences))
        p.map(RunIUPredAndUploadResults,ProteinSequences)
print('\n\nAnalysis Complete.\nTime Taken: %s' %(datetime.datetime.now()-start_time))
