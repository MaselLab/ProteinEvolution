import os, sys, json, csv, copy, mysql.connector,datetime

'''
Author : Sara Willis
Date   : May 17, 2019

The purpose of this script is to clear flagged Pfams from our MySQL database. When clearing Pfams from our dataset, there are a few steps that need to be followed to get the database in-line with the changes: 

   1) The pfams need to be removed from the pfam data table: PfamUIDs_EnsemblAndNCBI
   2) Any full genes that contain one of the removed pfams need to also be removed
   3) The pfams in the domain metrics table associated with the removed genes need to be removed
   4) Any scrambled sequences associated the removed genes should also be removed
   5) Species sets should be reassigned to pfams
   6) Pfams should be redated following species reassignment 
   7) Dating of full proteins should be redone
   8) Full genes should be reassigned full gene families. 
   9) Removed pfams should be pulled from the Alignment table

This script performs steps 2-4. The other steps can be performed using other scripts stored on Github

------------ WARNING ------------
This script deletes entries from preexisting tables! Make sure you have a recent backup made before running this script since you may want to return to your original dataset at a some point and no one wants to have to regenerate these tables from scratch!

As of May 17, 2019 Backups with the suffix 05062019 and 05072019 are tables with all pfam entries prior to cleanup

'''

##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################

# These are the databases where the protein, coding, protein metrics, coding metrics, and domain metrics
# are stored. These should all be related by the protein table UID and should come from the same source. i.e.
#     1) Multicellular
#     2) Ensembl
#     3) NCBI
#     4) Protists

# Data tables where sequences and metrics are stored
ProteinDatabase = 'Genomes_Multicellular_Protein'
CodingDatabase = 'Genomes_Multicellular_Coding'
ProteinMetricsDatabase = 'Genomes_Multicellular_ProteinMetrics'
DomainMetricsDatabase = 'Genomes_Multicellular_DomainMetrics'

# Ensembl and protist tables don't have Scrambled tables associated with them. For these, just set ScrambledDatabase = False
ScrambledDatabase = 'Genomes_Multicellular_ScrambledSequences'


# User-Specific MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# File with deleted pfams -- This should be a CSV file with a Pfam Accession to delete as the first entry per line
DeletedPfamsFile = 'PaulsListOfPfamsToRemove.csv'




##########################################################################################
#                                Program Executes Below                                  #
##########################################################################################

start_time = datetime.datetime.now()
print('\nProgram Executing\nCurrent Time: %s\n'%start_time)
sys.stdout.flush()

# A connection to the SQL database is established
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# The accessions for all pfams to delete are collected from the specified CSV file and
# stored in a list
DeletedPfams = []

with open(DeletedPfamsFile, 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for row in reader:
        # Sometimes Pfam accessions come version numbers which are attached following a period
        # since the pfams in our dataset don't have versions, we removed them if they exist
        # so we can accurately identify which removed pfams are in our dataset
        PfamUID = row[0].split('.')[0]
        DeletedPfams.append(PfamUID)


print('Extracting Data From Protein Database\n')
sys.stdout.flush()

# To delete the relevant information, we need the pfam UIDs and the protein table UID from the
# protein table
ExtractionStatement = "SELECT UID,PfamUID FROM "+ ProteinDatabase
mycursor.execute(ExtractionStatement)

print('Data Successfully Extracted\nTime Taken: %s\n\nParsing Results'%(datetime.datetime.now()-start_time))
sys.stdout.flush()

results = mycursor.fetchall()

n = 0
intermediate_time = datetime.datetime.now()
for result in results:
    ProteinTableUID, PfamUIDs = result
    
    # The pfam UIDs are stored as comma-delimited strings, so we partition them into a list for easy parsing
    PfamUIDs = PfamUIDs.split(',')
    for PfamUID in PfamUIDs:
        
        # If any of the deleted pfams show up in a protein, then that protein is deleted from the dataset
        if PfamUID in DeletedPfams:
            n += 1
            # The gene is deleted from the coding table
            DeleteFromCodingStatement = "DELETE FROM "+CodingDatabase+" WHERE UID = "+str(ProteinTableUID)
            mycursor.execute(DeleteFromCodingStatement)

            # The protein metrics table
            DeleteFromProteinMetricsStatement = "DELETE FROM "+ProteinMetricsDatabase+" WHERE ProteinTableUID="+str(ProteinTableUID)
            mycursor.execute(DeleteFromProteinMetricsStatement)

            # the protein table
            DeleteFromProteinStatement = "DELETE FROM "+ProteinDatabase + " WHERE UID="+str(ProteinTableUID)
            mycursor.execute(DeleteFromProteinStatement)

            # and any domains contained in that protein are deleted from the domain metrics table
            DeleteFromDomainMetricsStatement = "DELETE FROM "+DomainMetricsDatabase+" WHERE ProteinTableUID="+str(ProteinTableUID)
            mycursor.execute(DeleteFromDomainMetricsStatement)

            # and the scrambled versions of that protein are deleted from the scrambled table if they exist
            if ScrambledDatabase != False:
                DeleteFromScrambledStatement = "DELETE FROM " + ScrambledDatabase + " WHERE ProteinTableUID="+str(ProteinTableUID)
                mycursor.execute(DeleteFromScrambledStatement)

            cnx.commit()

            # The script notifies the user of the progress
            if n%1000 == 0:
                print('%s Pfam Results Removed\nTime Taken Since Last Update: %s'%(n,datetime.datetime.now()-intermediate_time))
                sys.stdout.flush()
                intermediate_time = datetime.datetime.now()
            # Since we found a pfam that qualifies for the deletion of the protein and all associated entries, we don't need to search
            # through that proteins pfams anymore since the entry was already removed, so we break the loop that was searching the pfams
            break
# Once the script ends, it notifies the user of how long it took and closes the sql connection
cnx.close()

print('Tables Successfully Cleaned\nTime Taken: %s\nPfam Entries Removed: %s'%(datetime.datetime.now()-start_time, n))
sys.stdout.flush()

        
