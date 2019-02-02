import json, sys, os, csv, mysql.connector, datetime
import numpy as np
from multiprocessing import Pool,Process


'''
The purpose of this script is to calculate the normalized indices of dispersion for protein sequences stored in a MySQL database. 


The functions used to calculate the index of dispersion were all created by Jason Bertram. The remainder of this script was written by Sara Willis.




IMPORTANT
---------
A note about speed: If you're updating a table and the protein sequence UID is not the primary key, make sure to convert the protein sequence UID to an index if the table is large. If the protein sequence UID is not indexed, the upload process could be *exceptionally* lengthy 

''' 


def LoadUserOptions():

    
    Database = '' # The name of the MySQL database where the relevant tables are stored
    user = '' # The username to access MySQL
    host = ''  # The IP address of the MySQL database
    password = '' # The password to access MySQL
    numberOfProcesses = 20 # The number of processes the user wishes to run in parallel
    hydroMap = 'FILVM' # The option in the hydrophobicity_map for which residues are counted as hydrophobic
    DispersionWindow = 6

    # Data Table Information
    ProteinDatabase = 'EnsemblGenomes_Protein_Complete' # The table where the protein sequences are stored
    ClusteringDatabase = 'EnsemblGenomes_ClusteringMetrics_Complete' # The table where the results should be saved -- can be the same or different from source table
    ProteinSeqTableUID_ColumnName = 'ProteinSeqTableUID' # The column name where protein UID should be stored in the destination table.
                                                         # If the destination table is the same as the source, this will be used to update the entries
    Trunc_Dispersion_FullGene_ColumnName = 'NormalizedIndexOfDispersion_Trunc_FILVM_FullGene' # Where to store the full gene dispersion index (trunc)
    Trunc_Dispersion_PfamOnly_ColumnName = 'NormalizedIndexOfDispersion_Trunc_FILVM_PfamOnly' # where to store the pfam dispersion index (trunc)
    AllPhases_Dispersion_FulLGene_ColumnName = 'NormalizedIndexOfDispersion_AllPhases_FILVM_FullGene' # Where to store the full gene dispersion index (all phases)
    AllPhases_Dispersion_PfamOnly_ColumnName = 'NormalizedIndexOfDispersion_AllPhases_FILVM_PfamOnly' # Were to store the pfam dispersion index (all phases)

    

    # Returns a dictionary so these options can be accessed later in the script
    return {'Database' : Database,
            'User' : user,
            'Host' : host,
            'Password' : password,
            'NumberOfProcesses' : numberOfProcesses,
            'ProteinDatabase' : ProteinDatabase,
            'ClusteringDatabase' : ClusteringDatabase,
            'DispersionWindow' : DispersionWindow,
            'HydroMap' : hydroMap,
            'ProteinSeqTableUID_ColumnName': ProteinSeqTableUID_ColumnName,
            'Trunc_Dispersion_FullGene_ColumnName': Trunc_Dispersion_FullGene_ColumnName,
            'Trunc_Dispersion_PfamOnly_ColumnName': Trunc_Dispersion_PfamOnly_ColumnName,
            'AllPhases_Dispersion_FullGene_ColumnName': AllPhases_Dispersion_FullGene_ColumnName,
            'AllPhases_Dispersion_PfamOnly_ColumnName': AllPhases_Dispersion_PfamOnly_ColumnName}




###################################################################################
########################## Jason's Clustering Functions ###########################
###################################################################################

'''
The following functions were written by Jason Bertram on Monday May 7 15:47:56 2018
All comments are direct quotes from his script.
'''
hydrophobicity_map={'FILVMW':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':1, 
                    'A':-1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':-1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1},
                    'FILVM':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':-1, 
                    'A':-1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':-1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1},
                    'AGWPY':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':0, 
                    'A':0,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':0,'H':-1,
                    'K':-1,'P':0,'S':-1,'T':-1,'Y':0,'X':0,'U':-1},
                    'FLIMVWAG':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':1, 
                    'A':1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1}}

#the original dispersion index with end truncation
#we set a window size of 6  (window=6) but this could in principle be messed around with
def fitness_dispersion_trunc(sequence_amino,window,hydro_map):
    length=float(len(sequence_amino))
    sequence_hydro=[hydrophobicity_map[hydro_map][_] for _ in sequence_amino]
    numblocks=np.floor(length/window)
    L=window*numblocks
    extraphases=length-L
    
    dispersion=0
    for phase in range(int(extraphases)+1):
        block_sums=np.array([np.sum(sequence_hydro[phase+window*_:phase+window*(_+1)]) for _ in range(int(numblocks))])
        M=np.sum(block_sums)
        if L**2==M**2 or numblocks<=1:
            dispersion=0
        else:
            dispersion=dispersion+(L-1)/((L**2-M**2)*(1-1/numblocks))*sum((block_sums-M/numblocks)**2)

    return dispersion/(1+extraphases)

#an alternative dispersion index which averages over all possible phases for the windowing    
def fitness_dispersion_all_phases(sequence_amino,window,hydro_map):
    length=float(len(sequence_amino))
    sequence_hydro=[hydrophobicity_map[hydro_map][_] for _ in sequence_amino]
        
    dispersion=0
    for phase in range(window):
        numblocks=np.floor((length-phase)/window)
        L=window*numblocks
        block_sums=np.array([np.sum(sequence_hydro[phase+window*_:phase+window*(_+1)]) for _ in range(int(numblocks))])
        M=np.sum(block_sums)
        if L**2==M**2 or numblocks<=1:
            dispersion=0
        else:
            dispersion=dispersion+(L-1)/((L**2-M**2)*(1-1/numblocks))*sum((block_sums-M/numblocks)**2)

    return dispersion/window


###################################################################################
###################################################################################
###################################################################################

'''
The following components of this script were written by Sara Willis on Monday January 7, 2019.

The script connects to the user's MySQL database, extracts protein sequences, runs them through Jason's clustering functions using the options specified at the beginning of this script and uploads the results to the user's database 

'''
# The program keeps track of the amount of time it takes for this script to execute. 
start_time = datetime.datetime.now()

# The program loads the user's connection information to access the MySQL database and connects
Database = LoadUserOptions()['Database']
User = LoadUserOptions()['User']
Host = LoadUserOptions()['Host']
Password = LoadUserOptions()['Password']
cnx = mysql.connector.connect(user = User,
                            password = Password,
                            host = Host,
                            database = Database)
# A cursor is defined so we can interact with the database
mycursor = cnx.cursor(buffered = True)



# This portion partitions the UIDs in the database so we can run clustering in parallel

# The maximum UID is extracted the from the data table so that the database can be divided into
# roughly equal-sized chunks
FindMaxUIDStatement = "SELECT MAX(UID) FROM %s"%LoadUserOptions()['ProteinDatabase']
mycursor.execute(FindMaxUIDStatement)
MaxUID = mycursor.fetchone()[0]

# A list is created to store the UIDs of the partitioned database. 
partitionedResults = []
sizeOfChunk =  MaxUID/LoadUserOptions()['NumberOfProcesses']
sizeOfChunk = int(sizeOfChunk)
chunks = list(range(0,MaxUID,sizeOfChunk))
if chunks[-1] < MaxUID:
    chunks.append(MaxUID)
for i in range(len(chunks)-1):
    start = chunks[i]
    stop = chunks[i+1]
    partitionedResults.append(list(range(start,stop)))
cnx.close()
batchNumber = 0
for element in partitionedResults:
    batchNumber += 1
    element.append(batchNumber)

# A function that opens a MySQL connection, extracts protein sequences, runs the sequences
# through Jason's clustering script, and uploads the results is defined. A function is created
# to perform this step so we can run it in parallel using the multiprocessing module.
def RunClustering(batch):

    # The script starts by assuming we're inserting entries into an empty table
    Update = False
    # If we're storing the results in the table where we extracted the protein sequences from, the program
    # is notified we need to use an update statement so we don't create new entries in the table
    if ProteinDatabase == ClusteringDatabase:
        Update = True
    # We establish a MySQL connection using the user options. This is in the function, opening a fresh connection
    # for each batch to avoid crosswiring.
    Database = LoadUserOptions()['Database']
    User = LoadUserOptions()['User']
    Host = LoadUserOptions()['Host']
    Password = LoadUserOptions()['Password']
    cnx = mysql.connector.connect(user = User,
                                  password = Password,
                                  host = Host,
                                  database = Database)
    mycursor = cnx.cursor(buffered = True)
    # The input for this function is a list with the following form:
    #     [UID_1, UID_2, ... , UID_N, BatchNumber]
    # The UIDs are used to extract entries from the datatable.
    # Batch Number is used to differentiate batches run in parallel allowing us
    # to save files from each process without overwriting files from other processes.
    
    batchNumber = batch[-1]
    UIDs = batch[:-1]

    for uid in UIDs:
        # We extract the UID, the protein sequence, and the pfam indices
        selectProteinSequenceStatement = "SELECT UID,ProteinSequence,PfamStart,PfamStop FROM %s WHERE UID = %s"%(LoadUserOptions()['ProteinDatabase'],uid)
        mycursor.execute(selectProteinSequenceStatement)
        result = mycursor.fetchone()
        # If there exists an entry with that UID, we continue. It's possible that there will be no entry with that
        # UID since entries may have been removed from the source data table. 
        if result != None:
            UID, ProteinSequence, Starts, Stops = result
            # The start/stop pfam indices are stored as lists.
            Starts = [int(i) for i in Starts.split(',')]
            Stops = [int(i) for i in Stops.split(',')]

            # We run both the trunc and allPhases function on our protein sequences using Jason's functions
            dispersion_trunc_fullGene = '%.4f'%fitness_dispersion_trunc(ProteinSequence,LoadUserOptions()['DispersionWindow'],LoadUserOptions()['HydroMap'])
            dispersion_allPhases_fullGene = '%.4f'%fitness_dispersion_all_phases(ProteinSequence,LoadUserOptions()['DispersionWindow'],LoadUserOptions()['HydroMap'])

            # If there are issues with a UID, we don't upload the results and instead save a file for later analysis
            problemUIDs = open('ProblemUIDs_%s.txt'%batchNumber,'w')
            UIDProblematic = False
            # We also want to calculate the dispersion indices for the Pfam domains, so we keep track of those values with lists
            PfamTrunc = []
            PfamAllPhases = []
            for i in range(len(Starts)):
                start = Starts[i]
                stop = Stops[i]
                try:
                    pfamSequence = ProteinSequence[start:stop]
                    dispersion_trunc_PfamOnly = '%.4f'%fitness_dispersion_trunc(pfamSequence, LoadUserOptions()['DispersionWindow'], LoadUserOptions()['HydroMap'])
                    dispersion_allPhases_PfamOnly = '%.4f'%fitness_dispersion_all_phases(pfamSequence,LoadUserOptions()['DispersionWindow'],LoadUserOptions()['HydroMap'])
                    PfamTrunc.append(dispersion_trunc_PfamOnly)
                    PfamAllPhases.append(dispersion_allPhases_PfamOnly)

                # It's possible there are some quality-control issues that need to be investigated. Specifically
                # Where either the protein sequence was truncated during the download process (it seems this
                # is particularly likely with downloads from the web-interface version of the Ensembl Biomart),
                # or there's some funkyness with the pfam annotations that make it so the indices appear outside
                # the body of the protein sequence. The UIDs of these problematic genes are stored in files so the
                # user knows when there was a problem and can go back and look at the issue.
                # The index error happens specifically when a Pfam index is called that is greater than the max
                # index of the protein sequence 
                except IndexError:
                    # if there was an error, the UID is saved to a file (the filename contains the batch number so files don't get
                    # overwritten if this function is being run in parallel). 
                    problemUIDs.write(UID)
                    # The function is then alerted there was an issue and the protein ceases to be processed. The results from
                    # these problematic genes are then not uploaded to the mySQL database. 
                    UIDProblematic = True
                    break
            # If no problems were found during the Pfam analysis, the results are uploaded to the destination table
            if UIDProblematic == False:
                # If the table we're uploading the results into isn't the same as the table we extracted the sequences from,
                # we need to check to see if the entries being processed already exist in the table or not
                if Update == False:
                    checkEntryExistence = "SELECT 1 FROM " + LoadUserOptions()['ClusteringDatabase'] +" WHERE "+ LoadUserOptions()['ProteinSeqTableUID_ColumnName'] +" = %s"
                    mycursor.execute(checkEntryExistence)
                    entry = mycursor.fetchone()
                    if entry == None:
                        
                        uploadResults = "INSERT INTO " + LoadUserOptions()['ClusteringDatabase'] + " (" ','.join([LoadUserOptions()['ProteinSeqTableUID_ColumnName'],LoadUserOptions()['Trunc_Dispersion_FullGene_ColumnName'],LoadUserOptions()['Trunc_Dispersion_PfamOnly_ColumnName'],LoadUserOptions()['AllPhases_Dispersion_FullGene_ColumnName'),LoadUserOptions()['AllPhases_Dispersion_PfamOnly_ColumnName']] +" VALUES) (%s,%s,%s,%s,%s)"
                    results = (UID,dispersion_trunc_fullGene,','.join(PfamTrunc),dispersion_allPhases_fullGene,','.join(PfamAllPhases))
                    mycursor.execute(uploadResults,results)
                    cnx.commit()
    problemUIDs.close()
    cnx.close()

# We process the full MySQL database in parallel with the number of processes defined by the user
if __name__ == '__main__':
    p = Pool(len(partitionedResults))
    p.map(RunClustering,partitionedResults)

# To finish, we consolidate all the problemUIDs into one file and delete the batch files
ProblemUIDsFile = open('ProblemUIDs.txt','w')
for i in range(batchNumber+1):
    if os.path.exists('ProblemUIDs_%s.txt'%i) == True:
        BatchProblemUIDsFile = open('ProblemUIDs_%s.txt'%i,'r')
        for line in BatchProblemUIDsFile:
            ProblemUIDsFile.write(line)
        os.remove(BatchProblemUIDsFile)
ProblemUIDsFile.close()  

print('Script Complete\nTime Taken: %s'%(datetime.datetime.now()-start_time))
