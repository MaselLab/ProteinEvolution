import os, sys, json, csv, mysql.connector, datetime
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool

'''
Author : Sara Willis
Date   : February 11, 2019
--------------------------

This script is used to run IUPred2 on proteins stored in a MySQL database.

**This is a slightly modified version of 'RunIUPred2.py'. This version ONLY deals with the raw IUPred output and doesn't bother with taking means

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
    IupredScript = '../iupred2a/iupred2a.py'                 # The path from this script's directory to the IUPred2 python script

    # User Options for IUPred
    # ------------------------
    DisorderType = 'long'                                   # Disorder options for IUPred. Can either be long, short, or glob
    CysteinesIncluded = True                                # True = include cysteines in input protein for analysis; False = excise cysteines

    # User's MySQL Database information
    # ------------------------
    Database = ''                                           # Name of user's database
    User = ''                                               # Username to access Fusion/MySQL
    Host = ''                                               # MySQL host
    Password = ''                                           # User's MySQL password
    DataTable = ''                                          # Name of the DataTable in the Database where sequences are stored
    ProteinColumnName = ''                                  # Name of the column where the user's proteins are stored in the data table
    UIDColumnName = ''                                      # Name of the column where the user's table uids are stored in the data table
    IupredRawScoreColumnName = ''                           # Name of the column where the user's raw IUPred scores are stored


    # The user-defined options are returned as a dictionary for use throughout the rest of the script
    return {'NumberOfProcesses': numberOfProcesses, 'PathToPython3': pathToPython3, 'IUPredScript': IupredScript, 'DisorderType': DisorderType, 'CysteinesIncluded' : CysteinesIncluded, 'RunAnchor' : RunAnchor, 'Database' : Database, 'User' : User, 'Host' : Host, 'Password' : Password, 'DataTable' : DataTable, 'ProteinColumnName' : ProteinColumnName, 'UIDColumnName' : UIDColumnName, 'IUPredRawScoreColumnName' : IupredRawScoreColumnName}



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

    # IUPred doesn't work if there's a stop codon symbol in the peptide sequence, so if there exists one, it's removed
    ProteinSequence = ProteinSequence.replace('*','')
    
    # IUPred reads in sequences from fasta files, so one is created for each protein sequence that is analyzed. It creates one with
    # the batch number attached because otherwise multiprocessing combined with tons of files that are all attempting to be written
    # and read with the same name is a no-go
    TempFastaFile = 'Temp_%s.fa'%batch
    IupredOutputFile = 'IUPredOutputFile_%s.txt'%batch

    # If a file with the desired name already exists, it's deleted
    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)

    CreateFastaFile(UID,ProteinSequence,TempFastaFile)

    if os.path.exists(IupredOutputFile) == True:
        os.remove(IupredOutputFile)

    os.system('%s %s %s %s >%s' %(UserOptions()['PathToPython3'],UserOptions()['IUPredScript'],TempFastaFile,UserOptions()['DisorderType'],IupredOutputFile))

    ISD_RawScore = []
    # Once Iupred is run, the output file is opened and parsed to return the relevant information to the user
    with open(IupredOutputFile,'r') as f:
        # Be default, the output from IUPred is a tab-delimited file
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            # If the row starts with the character #, it is human-readable and contains no data 
            if '#' in row or len(row) == 1:
                pass
            else:
                AANumber = int(row[0])
                AA = row[1] 
                ISDValue = row[2]
                ISD_RawScore.append(str(format(float(ISDValue),'.4f')))

    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)
    if os.path.exists(IupredOutputFile):
        os.remove(IupredOutputFile)

    return ','.join(ISD_RawScore)

    

# ----------------------------------------------------------------------------------
    

##### Partition Database #####

# The purpose of this function is to partition the database that is being processeed
# into roughly equal-sized chunks. The number of chunks is determined by the number
# of processes the user wishes to run

def PartitionDataset(MaxUID):
    
    partitionedUIDs = []

    sizeOfChunk = MaxUID/UserOptions()['NumberOfProcesses']
    sizeOfChunk = int(sizeOfChunk)

    chunks = list(range(0,MaxUID+1,sizeOfChunk))
    if chunks[-1] < MaxUID+1:
        chunks[-1] = MaxUID+1
    for i in range(len(chunks)-1):
        start = chunks[i]
        stop = chunks[i+1]
        partitionedUIDs.append(list(range(start,stop)))
    batchNumber = 0
    for element in partitionedUIDs:
        batchNumber += 1
        element.append(batchNumber)
    return partitionedUIDs

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
def RunIUPredAndUploadResults(ListOfUIDs):
    start_time = datetime.datetime.now()

    # The batch number is removed from the input list and stored
    batch = ListOfUIDs[-1]
    ListOfUIDs = ListOfUIDs[:-1]

    # A connection is then established with the MySQL server
    cnx = mysql.connector.connect(user = UserOptions()['User'],
                                  password = UserOptions()['Password'],
                                  host = UserOptions()['Host'],
                                  database = UserOptions()['Database'])
    # And a cursor is defined so the data tables can be interacted with 
    mycursor = cnx.cursor(buffered = True)
    # Each entry in the input list is then parsed and analyzed
    for uid in ListOfUIDs:

        selectProteinStatement = "SELECT " +UserOptions()['ProteinColumnName']+","+UserOptions()['UIDColumnName'] +" FROM "+UserOptions()['DataTable']+ " WHERE %s=%s"%(UserOptions()['UIDColumnName'],uid)
        mycursor.execute(selectProteinStatement)
        result = mycursor.fetchone()
        if result != None:
            protein= result[0]
            ISDRawScores = CalculateDisorder_IUPred2(uid, protein, batch)
            UploadStatement = "UPDATE " + UserOptions()["DataTable"] + " SET " + UserOptions()['IUPredRawScoreColumnName'] + "=%s WHERE " + UserOptions()['UIDColumnName'] + " = %s"
            UploadData = tuple((ISDRawScores,uid))
            mycursor.execute(UploadStatement,UploadData)
            cnx.commit()
    #print('Analysis Complete\nNumber of sequences processed: %s\nTime taken: %s'%(len(ListOfUIDs),(datetime.datetime.now()-start_time)))
    
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

PullMaxUID = "SELECT MAX(%s) FROM %s" %(UserOptions()['UIDColumnName'],UserOptions()['DataTable'])
mycursor.execute(PullMaxUID)
MaxUID = mycursor.fetchone()[0]
cnx.close()

PartitionedDatabase = PartitionDataset(MaxUID)

# Finally, we use multiprocessing to run the script in parallel as many
# times as the user specifies
if __name__ == '__main__':
    # In the event that the partitioning of the database results in a list that has fewer sublists
    # than the number of processes the user wants, we modify the pool number to match the length
    # of partitioned database to avoid any errors
    if len(PartitionedDatabase) == UserOptions()['NumberOfProcesses']:
        p = Pool(UserOptions()['NumberOfProcesses'])
        p.map(RunIUPredAndUploadResults,PartitionedDatabase)
    else:
        p = Pool(len(PartitionedDatabase))
        p.map(RunIUPredAndUploadResults,PartitionedDatabase)
print('\n\nAnalysis Complete.\nTime Taken: %s' %(datetime.datetime.now()-start_time))
