import os, json, sys, csv, mysql.connector, datetime, copy
import numpy as np

'''
Author : Sara Willis
Date   : February 19, 2019

The purpose of this script is to take proteins and their associated Pfam domains from a protein data table and to split those entries so there is a separate entry for each protein/pfam combination. Those entries are then uploaded into a new data table

The data table that the entries are being loaded into should already exist prior to running this script.

Note: The UIDs in the Protein database should be integers

'''

###########################################################################################
#                                User-Specific Information                                #        
###########################################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

ProteinDatabase = 'EnsemblGenomes_Protein_Complete'
Protein_UIDColumn = 'UID'
Protein_SpeciesUIDColumn = 'SpeciesUID'
Protein_PfamUIDColumn = 'PfamUID'
Protein_PfamStartColumn = 'PfamStart'
Protein_PfamStopColumn = 'PfamStop'

DomainMetricsDatabase = 'EnsemblGenomes_DomainMetrics_Complete_copy'
Domain_ProteinUIDColumn = 'ProteinTableUID'
Domain_SpeciesUIDColumn = 'SpeciesUID'
Domain_PfamUIDColumn = 'PfamUID'
Domain_LengthColumn = 'DomainLength'

###########################################################################################
#                             Program Executes Below                                      #         
###########################################################################################

# Establishes a connection to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
# A cursor allows us to interact with the database
mycursor = cnx.cursor(buffered = True)

# Extracts MaxUID from database so we can go through each entry, one at a time
ExtractionStatement = "SELECT Max(UID) FROM " + ProteinDatabase
mycursor.execute(ExtractionStatement)
MaxUID = mycursor.fetchone()[0]

# Loads the required columns so we can use them with ease
RequiredColumns_Protein = [Protein_UIDColumn,Protein_SpeciesUIDColumn,Protein_PfamUIDColumn,Protein_PfamStartColumn,Protein_PfamStopColumn]
RequiredColumns_Domain = [Domain_ProteinUIDColumn,Domain_SpeciesUIDColumn,Domain_PfamUIDColumn,Domain_LengthColumn]

# Goes through and extracts elements from the protein data table sequentially
for i in range(1, MaxUID+1):
    ExtractionStatement = "SELECT "+','.join(RequiredColumns_Protein)+" FROM "+ProteinDatabase+" WHERE "+Protein_UIDColumn+" = %s"%i
    mycursor.execute(ExtractionStatement)
    row = mycursor.fetchone()
    if row != None:
        UID,SpeciesUID,PfamUID,PfamStart,PfamStop = row
        PfamUID = PfamUID.split(',')
        PfamStart = [int(j) for j in PfamStart.split(',')]
        PfamStop = [int(j) for j in PfamStop.split(',')]
        # Then, for each pfam associated with a protein, an entry is made in the domain table for that protein/domain combo
        for j in range(0,len(PfamUID)):
            DomainLength = PfamStop[j] - PfamStart[j]
            UpdateStatement = "INSERT INTO "+DomainMetricsDatabase+' ('+ ','.join(RequiredColumns_Domain) + ") VALUES (%s,%s,%s,%s)"
            UpdateData = (UID,SpeciesUID,PfamUID[j],DomainLength)
            mycursor.execute(UpdateStatement,UpdateData)
            cnx.commit()


