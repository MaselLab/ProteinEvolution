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
The UIDs for the tables in question should be integers, otherwise this script will exit with an error

'''
#User-Specific MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

FullProteinDatabase = 'EnsemblGenomes_Protein_Complete'
DomainMetricsDatabase = 'EnsemblGenomes_DomainMetrics_Complete'




start_time = datetime.datetime.now()

print('\nBeginning Analysis\nCurrent Time: %s\n\n'%start_time)



cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)

mycursor = cnx.cursor(buffered = True)

ExtractionStatement = "SELECT MAX(UID) FROM "+FullProteinDatabase
mycursor.execute(ExtractionStatement)
MaxUID = mycursor.fetchone()[0]

n = 0
for i in range(1,MaxUID+1):
    SelectISDAndDomain = "Select UID,PfamUID,PfamStart,PfamStop,PfamMetricsTableUID,ISD_RawIUPredScores_CysIncluded FROM " +FullProteinDatabase + " WHERE UID = %s"%i
    mycursor.execute(SelectISDAndDomain)
    result = mycursor.fetchone()
    if result != None:
        UID, PfamUID,PfamStart, PfamStop, PfamMetricsTableUID, ISD = result
        if ISD == '':
            print(UID)
        else:
            ISD = [float(j) for j in ISD.split(',')]
            PfamStart = [int(j) for j in PfamStart.split(',')]
            PfamStop = [int(j) for j in PfamStop.split(',')]
            PfamMetricsTableUID = [int(j) for j in PfamMetricsTableUID.split(',')]

            for j in range(0,len(PfamStart)):
                try:
                    DomainISDList = ISD[PfamStart[j]:PfamStop[j]]
                    DomainMean = np.mean(DomainISDList)
                    DomainISD = '%.4f'%DomainMean
                    
                    updateStatement = "UPDATE "+ DomainMetricsDatabase +" SET MeanISD_IUPred2_WithCys=%s WHERE UID=%s"

                    updateData = (DomainISD,PfamMetricsTableUID[j])
                    mycursor.execute(updateStatement, updateData)
                    cnx.commit()
                except:
                    pass
cnx.close()

print('Analysis Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-start_time))
