import os, json, sys, csv, mysql.connector, datetime, copy, 
from scipy import stats
import numpy as np


'''
Author : Sara Willis
Date   : February 18, 2019

The purpose of this script is to pull coding sequences from a MySQL database and to check whether they pass a series of quality checks:

   1) The sequence starts with a start codon
   2) The sequence is a multiple of three in length
   3) There are no in-frame stop codons

Once it has been determined whether a sequence has passed these filters or not, the results are uploaded into a protein metrics table specified by the user
'''

# User information to access the MySQL database
Database = ''
User = ''
Host = ''
Password = ''

# This is where the coding sequences are stored 
CodingSequenceTable = 'EnsemblGenomes_Coding_Complete'

# This is where the protein sequence are stored. This is used to create a correspondance between the coding
# and protein table UIDs
ProteinSequenceTable = 'EnsemblGenomes_Protein_Complete'

# And this is the table where the results of the filtering are stored
ProteinMetricsTable = 'EnsemblGenomes_ProteinMetrics_Complete'


# When Verbose is set to True, the program keeps the user updated on the progress
# of the script
Verbose = True



##########################################################################
#                         Program Executes Below                         #
##########################################################################

# This defines a small function which takes in a coding sequence and checks to make sure the
# following three conditions are met:
#
#   1) The sequence starts with a start codon
#   2) The sequence is a multiple of three in length
#   3) There are no in-frame stop codons
#
# If any of those three conditions aren't met, the function returns False. Otherwise it returns
# True
def CodingSequenceQualityControl(CodingSequence):
    if CodingSequence[:3] != 'ATG':
        return False
    elif len(CodingSequence)%3 != 0:
        return False
    else:
        for i in range(0, len(CodingSequence)-3, 3):
            Codon = CodingSequence[i:i+3]
            if Codon == 'TAA' or Codon == 'TGA' or Codon =='TAG':
                return False
            else:
                return True


# The program keeps track of how long it takes to run this script
startTime = datetime.datetime.now()

# We establish a connection to MySQL
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
# And define a cursor so we can interact with the data tables
mycursor = cnx.cursor(buffered = True)


if Verbose == True:
    print('\n\nSelecting MaxUID\n\n')
    sys.stdout.flush()
# The max UID in the table is found so we can pull each sequence one by one from the table
mycursor.execute("SELECT MAX(UID) FROM " + CodingSequenceTable)
MaxUID = mycursor.fetchone()[0]

if Verbose == True:
    print('MaxUID Selected\nTime Taken: %s\n\nAssembling Coding Dictionary\n\n' %(datetime.datetime.now()-startTime))
    sys.stdout.flush()

# A dictionary where each coding sequence UID points to either True or False, depending on whether
# the sequence passed the quality filters
CodingDictionary = {}

for i in range(1,MaxUID+1):
    mycursor.execute("SELECT CodingSequence FROM "+CodingSequenceTable+" WHERE UID = %s",(i,))
    result = mycursor.fetchone()
    if result != None:
        CodingSequence = result[0]
        PassesQualityControl = CodingSequenceQualityControl(CodingSequence)
        CodingDictionary[i] = int(PassesQualityControl)

if Verbose == True:
    print('Coding Dictionary Assembled\nTime Elapsed: %s\n\nUpdating Tables\n\n'  %(datetime.datetime.now()-startTime))
    sys.stdout.flush()

# We then go through and pull the protein sequence UIDs and cross-check them with their coding sequences
mycursor.execute("SELECT MAX(UID) FROM "+ProteinSequenceTable)
MaxUID = mycursor.fetchone()[0]

for i in range(1,MaxUID+1):
    mycursor.execute("SELECT UID,CodingSeqTableUID FROM "+ProteinSequenceTable+" WHERE UID = %s",(i,))
    result = mycursor.fetchone()
    if result != None:
        # We use the coding sequence dictionary to pull out whether the sequence passed the filters or not
        # and update the protein metrics table. 
        UID,CodingSeqTableUID = result
        PassesQualityControl = CodingDictionary[CodingSeqTableUID]
        mycursor.execute("UPDATE "+ProteinMetricsTable+" SET PassedIUPredFilters=%s WHERE ProteinTableUID = %s",(PassesQualityControl,UID))
        cnx.commit()


print('Tables Updated\nTime Taken: %s' %(datetime.datetime.now()-startTime))
cnx.close()
        
