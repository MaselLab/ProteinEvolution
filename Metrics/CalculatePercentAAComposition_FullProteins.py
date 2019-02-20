import os, sys, json,csv, datetime,mysql.connector
import numpy as np

'''
Written : Friday January 11, 2019
Updated : Tuesday February 19, 2019
Author  : Sara Willis

This script is used to calculate the amino acid composition of protein sequences stored in a MySQL database. The script calculates the percent composition for each of the 20 amino acids, formats the values as a comma-delimited string, and updates a metrics table with the values. 

The metrics table should exist prior to running the script.

The UIDs associated with these tables should be integers, not strings

'''
############################################################################
#                       User-Specific Information                          #
############################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# Where the protein sequences are extracted from
ProteinDataTable = 'NCBIGenomes_Protein_Complete'
ProteinColumn = 'ProteinSequence'
UIDColumn = 'UID'

# Where the AA composition is saved
MetricsDataTable = 'NCBIGenomes_ProteinMetrics_Complete'
AACompositioncolumn = 'PercentAminoAcidComposition'
MetricsUIDColumn = 'ProteinTableUID'

Verbose = True

############################################################################
#                        Program Executes Below                            #
############################################################################

start_time = datetime.datetime.now()

if Verbose == True:
    print('Starting Analysis\nCurrent Time: %s'%start_time)
    sys.stdout.flush()

# Connects to the database 
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered=True)

# The way the amino acid percentages are saved must be uniform throughout the table, so we define
# a list of the amino acids so we count them in a systematic way
AminoAcids = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

mycursor.execute("SELECT MAX("+UIDColumn+") FROM %s"%ProteinDataTable)
MaxUID = mycursor.fetchone()[0]

# We iterate through each row in the data table one at a time
for i in range(1,MaxUID+1):
    SelectProteinsStatement = "SELECT "+','.join([UIDColumn,ProteinColumn])+" FROM "+ProteinDataTable+" WHERE "+UIDColumn+" = %s"
    UIDExtractionTuple = tuple((i,))
    mycursor.execute(SelectProteinsStatement,UIDExtractionTuple)
    Proteins = mycursor.fetchone()
    if Proteins != None:
        UID,ProteinSequence = Proteins
        ProteinLength = len(ProteinSequence)
        PercentAminoAcidList = []

        for AA in AminoAcids:
            RawNumber = ProteinSequence.count(AA)
            Percent = '%.5f'%(RawNumber/ProteinLength)
            PercentAminoAcidList.append(Percent)
            
        ReformattedPercent  = ','.join(PercentAminoAcidList)
        InsertionValues = (ReformattedPercent,UID)
        UpdateStatement = "UPDATE " + MetricsDataTable + " SET "+AACompositioncolumn+"=%s WHERE "+MetricsUIDColumn+"=%s"
        mycursor.execute(UpdateStatement,InsertionValues)
        cnx.commit()
    sys.exit(0)
