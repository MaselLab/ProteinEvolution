import csv, os, sys, json, mysql.connector, time, datetime, random
from multiprocessing import Pool, Process
from Bio import SeqIO
from Bio.Seq import Seq


'''
Author : Sara Willis
Date   : May 16, 2019

This script takes a set of protein sequences and generates a set of scrambled sequences by randomly sampling from the source protein without replacement. The sequences are then uploaded to a MySQL table to be stored along with the UID of the generating protein. The data table where the scrambled sequences will be stored should exist prior to running this script.
'''

##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################

# MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# Table with sequences to be scrambled
SourceTableName = 'Genomes_Multicellular_Protein'
ProteinColumn = 'ProteinSequence'
UIDColumn = 'UID'

# Destination table name
DestinationTableName = 'Genomes_Multicellular_ScrambledSequences'
ScrambledSequenceColumn = 'ScrambledSequence'
ProteinTableUIDColumn = 'ProteinTableUID'

# The number of scrambled sequences generated by each protein
NumberOfScrambledVersions = 1

# How "chatty" the program should be
Verbose = True


##########################################################################################
#                                Program Executes Below                                  #
##########################################################################################

start_time = datetime.datetime.now()

print('Program Executing\n\Generating %s Scrambled Sequences per Protein\nnCurrent Time: %s'%(NumberOfScrambledVersions,start_time))
sys.stdout.flush()

# Program establishes a connection to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# The max UID is extracted from the source data table so entries can be pulled one at a time
mycursor.execute("SELECT MAX("+UIDColumn+") FROM "+SourceTableName)
MaxUID = mycursor.fetchone()[0]

for uid in range(1,MaxUID):
    # The relevant entry is extracted from the database
    mycursor.execute("SELECT "+ProteinColumn+" FROM "+SourceTableName+" WHERE "+UIDColumn+"="+str(uid))
    ProteinSequence = mycursor.fetchone()
    if ProteinSequence != None:
        # If the user has selected Verbose, the program will notify the user each time a new 10000 entries have
        # been extracted
        if Verbose == True and uid%10000 ==0:
            print('%s Sequences Processed\nTime Elapsed: %s'%(uid,datetime.datetime.now()-start_time))
            sys.stdout.flush()
        # Stop codons are removed from the original sequence
        ProteinSequence = ProteinSequence[0].replace('*','')
        # And the number of scrambled sequences specified by the user are generated
        for i in range(NumberOfScrambledVersions):
            ScrambledSequence = ''.join(random.sample(ProteinSequence, len(ProteinSequence)))
            # Each scrambled sequence is then inserted into the relevant table
            mycursor.execute("INSERT INTO "+DestinationTableName+" ("+','.join([ProteinTableUIDColumn,ScrambledSequenceColumn])+") VALUES ("+','.join([str(uid),"'"+ScrambledSequence+"'"])+")")
            cnx.commit()
  

cnx.close()
print('Program Complete\nTime Taken: %s'%(datetime.datetime.now()-start_time))