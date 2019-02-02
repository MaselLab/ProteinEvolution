import os, json, sys, csv, mysql.connector, datetime
import numpy as np


'''
The purpose of this script is to calculate the GC content for all species in a MySQL database using their coding sequences. The results are then uploaded to be stored with the relevant species

'''

################################################################
#                   User-Specific Information                  #
################################################################

# MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''



################################################################
#                    Program Executes Below                    #
################################################################

start_time = datetime.datetime.now()
print('\nBeginning Analysis\nCurrent Time: %s\n\n'%start_time)

Databases = ['Ensembl','NCBI']

for SequenceDatabase in Databases:
    NucleotideSequenceDatabase = SequenceDatabase+"Genomes_Coding_Complete"

    # A connection is established to the SQL database
    cnx = mysql.connector.connect(user = User,
                                  password = Password,
                                  host = Host,
                                  database = Database)
    mycursor = cnx.cursor(buffered = True)

    # We want to iterate over all species that exist in the table, so we choose the max and
    # min species UID and iterate over that range
    mycursor.execute("SELECT MAX(SpeciesUID) FROM "+NucleotideSequenceDatabase)
    MaxSpeciesUID = mycursor.fetchone()[0]
    mycursor.execute("SELECT MIN(SpeciesUID) FROM "+NucleotideSequenceDatabase)
    MinSpeciesUID = mycursor.fetchone()[0]

    for i in range(MinSpeciesUID,MaxSpeciesUID):
        # GC content is the fraction of G's and C's to all nucleotides
        FullLength = 0
        SumOfGsAndCs = 0
        
        mycursor.execute("SELECT CodingSequence FROM " + NucleotideSequenceDatabase + " WHERE SpeciesUID = %s"%i)
        results = mycursor.fetchall()
        if results == []:
            pass
        else:
            for result in results:
                NucleotideSequence = result[0]
                FullLength += len(NucleotideSequence)
                SumOfGsAndCs += (NucleotideSequence.count('C') + NucleotideSequence.count('G'))

            FractionOfGC = '%.4f'%(SumOfGsAndCs/FullLength)
            # Once the GC content is found, we store it in a species database for later use
            mycursor.execute("UPDATE SpeciesList SET GCContent=%s WHERE SpeciesUID = %s",(FractionOfGC,i))
            cnx.commit()
cnx.close()
    

