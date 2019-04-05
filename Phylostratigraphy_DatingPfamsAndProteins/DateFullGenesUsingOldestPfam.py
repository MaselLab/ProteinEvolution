import os, sys, json, mysql.connector, datetime, csv, copy
import numpy as np


'''
Author       : Sara Willis
Date         : February 1, 2019
Updated      : April 5, 2019
Update Notes : Changed so MySQL column names defined only once allowing for greater ease-of-use

---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- 

This script is used to date full proteins by assigning the age of the oldest Pfam that it contains. 

In order to run this script, the following databases are required:

   - Two databases containing protein sequences and the pfams that occur within them (one for NCBI, one for Ensembl)
   - Two databases where the ages of the protein sequences will be stores (one for NCBI, one for Ensembl)
   - A database where pfams are stored with their ages. The pfam UID should be the primary key in this database

The following databases are what have been used:

In the database PFAMphylostratigraphy, there are four data tables that contain full proteins and the metrics associated with them:

   1) EnsemblGenomes_Proteins_Complete
   2) EnsemblGenomes_ProteinMetrics_Complete
   3) NCBIGenomes_Proteins_Complete
   4) NCBIGenomes_ProteinMetrics_Complete

The Proteins data table contains a variety of data that includes the UIDs of the pfams that show up in the body of the protein. The ProteinMetrics data tables contain the averages of the metrics that have been calculated for each protein as well as the age of the oldest Pfam that occurs in the protein. The UID of the ProteinMetrics data table is shared with the Proteins data table. 

For each of the source databases (Ensemble and NCBI), the PfamUIDs are extracted from the Proteins database for each protein, the oldest is found, and the age of the oldest is used to update the ProteinMetrics data table. 

In PFAMphylostratigraphy, the data table PfamUIDsTable_EnsemblAndNCBI contains all PfamUIDs and the ages that have been assigned to them.
'''


########################################################################
########                    User-Supplied Data                  ########
########################################################################

# MySQL Connection Information 
Database = ''
User = ''
Host = ''
Password = ''

# Protein Data Tables (contains relevant Pfam UIDs)
EnsemblProteinDataTable = 'EnsemblGenomes_Protein_Complete'
EnsemblUIDColumn = 'UID'
EnsemblPfamUIDColumn = 'PfamUID'

NCBIProteinDataTable = 'NCBIGenomes_Protein_Complete'
NCBIUIDColumn = 'UID'
NCBIPfamUIDColumn = 'PfamUID'

# Protein Metrics Data Table (where the protein ages should be stored)
EnsemblProteinMetricsDataTable = 'EnsemblGenomes_ProteinMetrics_Complete'
EnsemblAgeOfOldestPfamColumn = 'AgeOfOldestPfam_LUCA_LECA_Updated'
EnsemblProteinTableUID = 'ProteinTableUID'

NCBIProteinMetricsDataTable = 'NCBIGenomes_ProteinMetrics_Complete'
NCBIAgeOfOldestPfamColumn = 'AgeOfOldestPfam_LUCA_LECA_Updated'
NCBIProteinTableUID = 'ProteinTableUID'

# Pfam Data Table
PfamDataTable = 'PfamUIDsTable_EnsemblAndNCBI'
PfamUIDColumn = 'PfamUID'
PfamAgeColumn = 'Age_Oldest_MY'

Verbose = True


########################################################################
########                  Program Executes Below                ########
########################################################################


# First the program establishes a connection to the SQL database
startTime = datetime.datetime.now()
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)


# The protein table UID and the pfam UIDs associated with that protein are extracted from the Ensembl database
if Verbose == True:
    print('Extracting Entries from Ensembl Protein Data Table\n\n')
    sys.stdout.flush()

SelectProteins_Ensembl_Statement = "SELECT "+','.join([EnsemblUIDColumn,EnsemblPfamUIDColumn])+" FROM " + EnsemblProteinDataTable
mycursor.execute(SelectProteins_Ensembl_Statement)
EnsemblResults = mycursor.fetchall()

# The process of dating the genes starts once the extraction is complete
if Verbose == True:
    print('Ensembl Entries Extracted\nTime Taken: %s\n\nDating Ensembl Proteins\n\n'%(datetime.datetime.now()-startTime))
    sys.stdout.flush()
currentTime = datetime.datetime.now()
for result in EnsemblResults:
    UID = result[0]
    Pfams = result[1].split(',')
    OldestPfamAge = 0
    # The pfam ages associated with that protein are then extracted from the pfam table.
    # The ages are iterated through and the oldest is selected
    for pfam in Pfams:
        SelectPfamStatement = "SELECT "+PfamAgeColumn+" FROM " + PfamDataTable + " WHERE "+PfamUIDColumn+" = "+pfam
        mycursor.execute(SelectPfamStatement)
        Age = mycursor.fetchone()[0]
        if Age >= OldestPfamAge:
            OldestPfamAge = copy.deepcopy(Age)
    # Once the oldest pfam is found, that age is assigned to the protein. That age is then uploaded
    # to the protein metrics table
    UpdateProteinFullAgeStatement = "UPDATE " + EnsemblProteinMetricsDataTable + " SET "+EnsemblAgeOfOldestPfamColumn+"="+OldestPfamAge+" WHERE "+EnsemblProteinTableUID+" ="+str(UID)
    mycursor.execute(UpdateProteinFullAgeStatement)
    cnx.commit()

if Verbose == True:
    print('Ensembl Protein Table Updated\nTime Taken: %s\n\nExtracting Entries from NCBI Protein Data Table\n\n'%(datetime.datetime.now()-currentTime))
    sys.stdout.flush()
currentTime = datetime.datetime.now()

# Once the Ensembl proteins are dated, the NCBI proteins go through the same process
SelectProteins_NCBI_Statement = "SELECT "+','.join([NCBIUIDColumn,NCBIPfamUIDColumn])+" FROM " + NCBIProteinDataTable
mycursor.execute(SelectProteins_NCBI_Statement)
NCBIResults = mycursor.fetchall()
print('NCBI Entries Extracted\nTime Taken: %s\n\nDating NCBI Proteins' %(datetime.datetime.now()-currentTime))
currentTime = datetime.datetime.now()
for result in NCBIResults:
    UID = result[0]
    Pfams = result[1].split(',')
    OldestPfamAge = 0
    for pfam in Pfams:
        SelectPfamStatement = "SELECT "+PfamAgeColumn+" FROM " + PfamDataTable + " WHERE "+PfamUIDColumn+" = "+pfam
        mycursor.execute(SelectPfamStatement)
        Age = mycursor.fetchone()[0]
        if Age >= OldestPfamAge:
            OldestPfamAge = copy.deepcopy(Age)
    UpdateProteinFullAgeStatement = "UPDATE " + NCBIProteinMetricsDataTable + " SET "+NCBIAgeOfOldestPfamColumn+"="+OldestPfamAge+" WHERE "+NCBIProteinTableUID+" = "+str(UID)
    mycursor.execute(UpdateProteinFullAgeStatement)
    cnx.commit()
if Verbose == True:
    print('NCBI Protein Table Updated\nTime Take: %s\n\nProgram Complete!\n\nTotal Time Taken: %s'%(datetime.datetime.now()-currentTime, datetime.datetime.now()-startTime))
    sys.stdout.flush()

cnx.close()

print('\nScript Complete\nTime Taken: %s'%(datetime.datetime.now()-startTime))
