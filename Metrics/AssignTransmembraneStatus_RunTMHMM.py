

import tempfile
import os, sys, json, csv, mysql.connector, datetime
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool

import os, sys, json, csv, re, subprocess, mysql.connector, datetime


'''

Authors : Sara Willis and Jennifer James
Date    : May 31, 2019

This script is used to run tmhmm on protein sequences stored in a MySQL database, and upload the output.
TMHMM is designed for the prediction of transmembrane helices in proteins, and is designed to be run on full protein sequences. 
More details at http://www.cbs.dtu.dk/services/TMHMM/

This script assumes the columns required for the tmhm output already exits in the table being updated (see tmhmm_run for column names)


Dependencies:
--------------
--> anaconda3 (specifically python3) : https://www.anaconda.com/download/
--> mysql.connector : https://pynative.com/install-mysql-connector-python/
--> Biopython : https://anaconda.org/conda-forge/biopython
--> tmhmm : installed locally
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
    numberOfProcesses = 10

    # User Paths
    # ------------------------
    pathToPython3 = '~/anaconda3/bin/python3'               # This is the path to the python3 executable. IUPred2 requires Python3 to run
    tmhmm_path = "~/tmhmm-2.0c/bin/tmhmm"                   # full user path + extension for tmhmm
    homedir = ''                              # user specified location of tempfile and tmhmm output 

    # User's MySQL Database information
    # ------------------------
    Database = ''                      # Name of user's database
    User = ''                                               # Username to access Fusion/MySQL
    Host = ''                                      # MySQL host
    Password = ''                                           # User's MySQL password
    DataTable = 'NCBIGenomes_Protein_Complete'              # Name of the DataTable in the Database where sequences are stored
    ResultsTable = 'NCBIGenomes_ProteinMetrics_Complete'    # Name of the DataTable in the Database where the results are stored
    ProteinColumnName = 'ProteinSequence'                   # Name of the column where the user's proteins are stored in the data table
    UIDColumnName = 'UID'                                   # Name of the column where the user's table uids are stored in the data table

    # The user-defined options are returned as a dictionary for use throughout the rest of the script
    return {'NumberOfProcesses': numberOfProcesses, 'PathToPython3': pathToPython3, 'tmhmm_path': tmhmm_path, 'homedir': homedir, 'Database' : Database, 'User' : User, 'Host' : Host, 'Password' : Password, 'DataTable' : DataTable, 'ProteinColumnName' : ProteinColumnName, 'UIDColumnName' : UIDColumnName, 'ResultsTable' : ResultsTable}


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

### Run tmhmm ###

'''
The purpose of this submodule is to call tmhmm, and pipe the output back to us so we can update the appropriate column in our results table.

TMHMM output overview:
In the short output format one line is produced for each protein with no
graphics. Each line starts with the sequence identifier and then these
fields:

"ExpAA=": Exp number of AAs in TMHs: The expected number of amino acids intransmembrane
helices. If this number is larger than 18 it is very likely to be a transmembrane
protein (OR have a signal peptide)

"First60=": Exp number, first 60 AAs: The expected number of amino acids in transmembrane
helices in the first 60 amino acids of the protein. If this number more
than a few, you should be warned that a predicted transmembrane helix in
the N-term could be a signal peptide.

"PredHel=": The number of predicted transmembrane helices by N-best, the overall most probable structure.

"Topology=": The topology predicted by N-best.
Example output: i7-29o44-66i87-109o
The topology is given as the position of the transmembrane helices separated
by 'i' if the loop is on the inside or 'o' if it is on the outside. The
above example '<tt>i7-29o44-66i87-109o' means that it starts on the inside,
has a predicted TMH at position 7 to 29, the outside, then a TMH at position
44-66 etc.

'''

def tmhmm_run(ProteinSequences):

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

	for ProteinSequence in ProteinSequences:
		### first clear fasta file
		fp = tempfile.NamedTemporaryFile(dir = UserOptions()['homedir'], delete = False)
		
		### fasta file with header as UID
		print(fp.name)
		CreateFastaFile(ProteinSequence[0], ProteinSequence[1], fp.name)
		### call to subprocess to run tmhmm on shell
		tmhmm_results = subprocess.Popen(UserOptions()['tmhmm_path']+" "+fp.name+" -short", stdout=subprocess.PIPE, shell = True)
		### get output	
		out, err = tmhmm_results.communicate()
		total_tmhmm_results =  str(out)[2:-3] ### cleaned up tmhmm results, single line tab separated
		ExpAA = total_tmhmm_results.split('\\t')[2].split('=')[1]
		First60 = total_tmhmm_results.split('\\t')[3].split('=')[1]
		PredHel = total_tmhmm_results.split('\\t')[4].split('=')[1]
		Topology = total_tmhmm_results.split('\\t')[5].split('=')[1]
		Topology = Topology.strip('\\n')

		update_statement = "INSERT INTO "+UserOptions()['ResultsTable']+" SET ExpAA = "+str(ExpAA)+", First60 = "+str(First60)+", PredHel = "+str(PredHel)+", Topology = '" +str(Topology) + "' WHERE ProteinTableUID = "+str(ProteinSequence[0])

		print(update_statement)

		mycursor.execute(update_statement)
		cnx.commit()
		os.remove(fp)
	cnx.close()


# ----------------------------------------------------------------------------------


####################################################################################
############################# Program Executes Below ###############################
####################################################################################


start_time = datetime.datetime.now()

print('Beginning tmhmm calculation\nNumber of Processes Being Used: %s\nCurrent Time: %s'%(UserOptions()['NumberOfProcesses'],start_time))

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
### tuples of UID and Protein sequence
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
        p.map(tmhmm_run,ProteinSequences)
    else:
        p = Pool(len(ProteinSequences))
        p.map(tmhmm_run,ProteinSequences)
print('\n\nAnalysis Complete.\nTime Taken: %s' %(datetime.datetime.now()-start_time))

