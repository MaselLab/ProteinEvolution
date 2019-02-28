import os, sys, csv, mysql.connector
from FindMitochondrialGenes_Submodules import LoadConnectionInformation

''''
Author : Sara Willis
Date   : Wednesday February 27, 2019



The purpose of this file is to remove identified mitochondrial genes from the Ensembl MySQL dataset. 

The transcript IDs were initially pulled from Ensembl's FTP cDNA files for the species in the database.

The coding sequences from the Ensembl database are pulled and the transcript IDs are compared with those identified as mitochondrial. The entries' UIDs are then used to upload the mitochondrial transcripts to a backup table and the entries are deleted from the main datatable.
'''




# A dictionary is created to store the datatable entries
GenesDatabase = {}

# An empty list is defined to store the UIDs of the sequences identified as mitochondrial
MitochondrialGenesToRemove = []

# The name of the flat file containing the transcript IDs corresponding to mitochondrial sequences is defined
mitochondrialGenesFile = 'EnsemblChloroplastGenes.txt'

# The information to log into the MySQL database is then loaded. It is stored in a separate file to protect
# the user's login information
Database = ''
User = ''
Host = ''
Password = ''

# We then connect to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)

# And create a cursor so we can interact with the table. 
mycursor = cnx.cursor(buffered = True)

# We define a statement to pull the UIDs and Transcript IDs from the coding table. We start
# with the coding table because it is smaller than the protein table and so is less computationally expensive
# to work with
pullGenesStatement = "SELECT UID,TranscriptID FROM EnsemblGenomes_Protein_Complete"

# and we pull the entries from the table
mycursor.execute(pullGenesStatement)

# Each entry pulled is sorted through
for entry in mycursor:
    UID = entry[0]
    TranscriptID = entry[1].split('.')[0]

    # We then add the sequences to the defined dictionary.
    # There shouldn't be duplicate TranscriptIDs, but to account for if there are
    # we make each entry in the dictionary a list of UIDs associated with the transcript UID so if more than one
    # shows up, we append it to the list so we don't neglect it
    if TranscriptID not in GenesDatabase:
        GenesDatabase[TranscriptID] = [UID]
    else:
        GenesDatabase[TranscriptID].append(UID)

        
# We then open the mitochondrial genes file to extract the transcript IDs
with open(mitochondrialGenesFile, 'r') as f:
    # We define a reader
    reader = csv.reader(f)
    # And go through the file row by row
    for row in reader:
        # The transcript ID is the only entry in each row. 
        transcriptID = row[0]
        # We then search the transcript ID against our dictionary that we created for our database.
        if transcriptID in GenesDatabase:
            # The UIDs associated with the transcriptID are added to the list of UIDs to remove
            MitochondrialGenesToRemove += GenesDatabase[transcriptID]

# To ensure the UIDs in the list are only represented once, we convert the list to a set and back
MitochondrialGenesToRemove = list(set(MitochondrialGenesToRemove))

# We then define two empty dictionaries to store the rows associated with the entries being removed. One for the coding
# sequences, one for the protein sequences
ToUpload_Coding = {}
ToUpload_Protein = {}

# Then the UIDs that will be removed are used to pull rows from the two datatables.
for UID in MitochondrialGenesToRemove:
    
    # The extraction statement is defined for the coding table
    pullGenesWithUIDStatement = "SELECT * FROM EnsemblGenomes_Coding_Complete WHERE UID =%s"%UID
    # The extraction statement is executed
    mycursor.execute(pullGenesWithUIDStatement)
    # and the result is pulled using the fetchone statement
    result = mycursor.fetchone()
    # The row is then stored in the the Coding dictionary using its UID as the key
    ToUpload_Coding[result[0]] = result
    
    # The same process is done for the protein table then
    pullProteinsWithUIDStatement = "SELECT * FROM EnsemblGenomes_Protein_Complete WHERE UID=%s"%UID
    mycursor.execute(pullProteinsWithUIDStatement)
    result = mycursor.fetchone()
    ToUpload_Protein[result[0]] = result

# Once the dictionaries are complete, the entries need to be uploaded to the datatable. We start with the coding table
# Each key is searched

for key in ToUpload_Coding:
    # An insertion statement is defined
    uploadGenesStatement = "INSERT INTO RemovedMitochondrialGenes_Coding (UID,Repository,SpeciesUID,NewickSpecies,GeneID,TranscriptID,ProteinID,TranscriptSupportLevel,CodingSequence) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    # And the row is uploaded to the table
    mycursor.execute(uploadGenesStatement,ToUpload_Coding[key])

    # The sequences are then removed from the full database
    # The deletion statement is defined
    deleteCodingSequenceStatement = "DELETE FROM EnsemblGenomes_Coding_Complete WHERE UID = %s"%UID
    # And is executed
    mycursor.execute(deleteCodingSequenceStatement)

    # And the changes are then committed
    cnx.commit()

# The same process is executed for the protein table
for key in ToUpload_Protein:
    
    uploadProteinsStatement = "INSERT INTO RemovedMitochondrialGenes_Protein (UID,CodingSeqTableUID,SpeciesUID,NewickSpecies,GeneID,TranscriptID,ProteinID,ProteinSequence,PfamUID,PfamStart,PfamStop,InterProScanEval,NoCysMeanISD,NoCysISD_RawIUPredScores,PfamDomainMeanISDValues) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    try:
        mycursor.execute(uploadProteinsStatement,ToUpload_Protein[key])
        deleteProteinSequenceStatement = "DELETE FROM EnsemblGenomes_Protein_Complete WHERE CodingSeqTableUID=%s"%UID
        mycursor.execute(deleteProteinSequenceStatement)
        cnx.commit()
    except:
        try:
            deleteProteinSequenceStatement = "DELETE FROM EnsemblGenomes_Protein_Complete WHERE CodingSeqTableUID=%s"%UID
            mycursor.execute(deleteProteinSequenceStatement)
            cnx.commit()
        except:
            pass

# The mysql connection is then terminated before the program ends
cnx.close()
