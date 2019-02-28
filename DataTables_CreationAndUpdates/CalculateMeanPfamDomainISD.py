import os, json, sys, csv, mysql.connector, datetime
import numpy as np


'''
Author : Sara Willis
Date   : February 11, 2019
--------------------------


This script is intended to pull preexisiting raw IUpred scores from a full-protein run from a MySQL data table and to find the mean disorder prediction for each of its Pfam domains. The results get uploaded to a domain data table

To use this script, the IUPred scores should be stored in a database with comma-delimited strings of pfams IDs, their starting and stopping indices, and a list of UIDs associated with those pfam domains in a domain data table

NOTE: 
-----
The UIDs for the tables in question should be integers

'''

################################################################################
#                        User-Specific Information                             #
################################################################################

# MySQL connection information
Database = ''
User = ''
Host = ''
Password = ''

# Database where raw IUPred scores are stored
FullProteinDatabase = 'EnsemblGenomes_Protein_Complete'
# Columns in Protein Database 
ProteinUIDColumn = 'UID'
PfamUIDColumn = 'PfamUID'
PfamStartColumn = 'PfamStart'
PfamStopColumn = 'PfamStop'
PfamMetricsTableUIDColumn = 'PfamMetricsTableUID'
RawIUPredScoresColumn = 'ISD_RawIUPredScores_CysIncluded'

# Database where the pfam domain metrics are stored
DomainMetricsDatabase = 'EnsemblGenomes_DomainMetrics_Complete'
#Columns in Domain Database
DomainMeanISDColumn = 'MeanISD_IUPred2_WithCys'
DomainTableUID = 'UID'

# When Verbose is set to True, the program will keep the user informed on the script's progress
Verbose = True


################################################################################
#                           Program Executes Below                             #
################################################################################


start_time = datetime.datetime.now()

if Verbose == True:
    print('\nBeginning Analysis\nCurrent Time: %s\n\n'%start_time)
    sys.stdout.flush()
    
# Connects to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We start by extracting the max UID from the protein database so we can iterate through each row in the database
ExtractionStatement = "SELECT MAX("+ProteinUIDColumn+") FROM "+FullProteinDatabase
mycursor.execute(ExtractionStatement)
MaxUID = mycursor.fetchone()[0]

# We're going to keep an eye on whether any entries in the data table didn't have an ISD score associated with it. This will be reported to the user
# at the end of the script to keep the user aware of the status of their data table. This will only be reported if Verbose is set to True
AnISDValueDoesNotExist = False

# We also want to keep track of whether there is any funkyness with our Pfam Domains. If there is, we keep note of it. This can happen if there's
# bad Pfam annotation or the full protein sequence was not stored correctly in the protein database. This is particularly common if a sequence
# is very long and was manually uploaded via navicat. When this happens, it's possible that a pfam index will point to a section of the protein
# sequence that was not uploaded, i.e. an index will be greater than the length of the present peptide sequence. We'll keep track of the problematic
# UIDs and, if the user has set Verbose to True, the user will be returned the problematic UIDs
ProblematicPfamIndicesExist = False
ProblematicUIDs = set([])

# We iterate through the UIDs from 1 to MaxUID, extracting every row. The reason we add 1 to the range is because the last element in range(1,n)
# is n-1
for i in range(1,MaxUID+1):
    # We extract all information from the columns, whose names are supplied in the user-specific portion of this script
    SelectISDAndDomain = "Select "+','.join([ProteinUIDColumn,PfamUIDColumn,PfamStartColumn,PfamStopColumn,PfamMetricsTableUIDColumn,RawIUPredScoresColumn])+" FROM " +FullProteinDatabase + " WHERE "+ProteinUIDColumn+" = %s"
    UIDTuple = tuple((i,))
    mycursor.execute(SelectISDAndDomain,UIDTuple)
    result = mycursor.fetchone()
    # If entries are removed from the data table, the entry associated with that UID no longer exists. We avoid these by making sure we
    # have extracted nothing from the table before we proceed
    if result != None:
        UID, PfamUID,PfamStart, PfamStop, PfamMetricsTableUID, ISD = result
        if ISD == '':
            AnISDValueDoesNotExist = True
        else:
            # We split up pfams, their UIDs, and their indices within the body of the protein into lists. These values are used
            # to compute the mean ISD for the domains
            ISD = [float(j) for j in ISD.split(',')]
            PfamStart = [int(j) for j in PfamStart.split(',')]
            PfamStop = [int(j) for j in PfamStop.split(',')]
            PfamMetricsTableUID = [int(j) for j in PfamMetricsTableUID.split(',')]

            # For each pfamUID, we extract the scores associated with its indices within the protein from the IUPred raw string
            # and average them to get the mean ISD for that particular domain
            for j in range(0,len(PfamStart)):
                try:
                    # Calculates domain ISD 
                    DomainISDList = ISD[PfamStart[j]:PfamStop[j]]
                    DomainMean = np.mean(DomainISDList)
                    DomainISD = '%.4f'%DomainMean

                    # and updates domain metrics table
                    updateStatement = "UPDATE "+ DomainMetricsDatabase +" SET "+DomainMeanISDColumn+"=%s WHERE "+DomainTableUID+"=%s"
                    updateData = (DomainISD,PfamMetricsTableUID[j])
                    mycursor.execute(updateStatement, updateData)
                    cnx.commit()
                    
                except:
                    ProblematicPfamIndicesExist = True
                    ProblematicUIDs.add(i)

cnx.close()

# If the user has selected Verbose = True, the user is given a rundown on the programs execution.
if Verbose == True:
    
    print('Analysis Complete\n\n')
    if AnISDValueDoesNotExist == True:
        print('Warning: A UID in your protein table does not have ISD data associated with it\n\n')
    if ProblematicPfamIndicesExist == True:
        print('Warning: There were index errors in Pfam analysis. Writing UIDs of problematic proteins to file\n')
        ProblemUIDsFilename = 'UIDsWithIndexErrors.txt'
        ProblemUIDs = open(ProblemUIDsFilename,'w')
        for i in ProblematicUIDs:
            ProblemUIDs.write('%s\n'%i)
        ProblemUIDs.close()
        print('File Written\nFilename: %s\n\n'%ProblemUIDsFilename)
    print('Total time taken: %s'%(datetime.datetime.now()-start_time))
