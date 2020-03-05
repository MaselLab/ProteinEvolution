import os, json, sys, csv, mysql.connector, datetime
import numpy as np


'''
Author : Jenny James and Sara Willis
Date   : January, 2020
--------------------------
This script will calculate the mean amino acid composition scores in % for pfam domains in the database, and upload to a preexisting database column.

NOTE to database users: this script uploads to a mysql column that will be renamed- an older 'PercentAminoAcidComposition' column was replaced as this run
                        did not complete correctly. 
                        
'''
#User-Specific MySQL Connection Information
Database = 'PFAMphylostratigraphy'
User = ''
Host = '127.0.0.1'
Password = ''

FullProteinDatabase = 'NCBIGenomes_Protein_Complete'
DomainMetricsDatabase = 'NCBIGenomes_DomainMetrics_Complete'
AACompositioncolumn = 'PercentAminoAcidCompositionCheck'


start_time = datetime.datetime.now()

print('\nBeginning Analysis\nCurrent Time: %s\n\n'%start_time)



cnx = mysql.connector.connect(user = User,
							  password = Password,
							  host = Host,
							  database = Database)

mycursor = cnx.cursor(buffered = True)

# The way the amino acid percentages are saved must be uniform throughout the table, so we define
# a list of the amino acids so we count them in a systematic way
AminoAcids = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

ExtractionStatement = "SELECT MAX(UID) FROM "+FullProteinDatabase
mycursor.execute(ExtractionStatement)
MaxUID = mycursor.fetchone()[0]

n = 0
for i in range(1,MaxUID+1):
	SelectISDAndDomain = "Select UID,PfamUID,PfamStart,PfamStop,PfamMetricsTableUID,ProteinSequence FROM " +FullProteinDatabase + " WHERE UID = %s"%i
	mycursor.execute(SelectISDAndDomain)
	result = mycursor.fetchone()
	if result != None:
		UID, PfamUID,PfamStart, PfamStop, PfamMetricsTableUID, Sequence = result
		if Sequence == '':
			print(UID)
		
		else:
			PfamStart = [int(j) for j in PfamStart.split(',')]
			PfamStop = [int(j) for j in PfamStop.split(',')]
			PfamMetricsTableUID = [int(j) for j in PfamMetricsTableUID.split(',')]

			for j in range(0,len(PfamStart)):
				try:
				
					### because positions use inclusive counting	
					ProteinLength = int(PfamStop[j])+1-int(PfamStart[j])
					PercentAminoAcidList = []
					
					### because strings are 0 starting, and slicing is inclusive to start and exclusive to stop
					DomainSequence = Sequence[PfamStart[j]-1:PfamStop[j]]
					
					for AA in AminoAcids:
					
						RawNumber = DomainSequence.count(AA)
						Percent = '%.5f'%(RawNumber/ProteinLength)
						PercentAminoAcidList.append(Percent)
				
					ReformattedPercent  = ','.join(PercentAminoAcidList)
					InsertionValues = (ReformattedPercent,PfamMetricsTableUID[j])   
					
					UpdateStatement = "UPDATE " + DomainMetricsDatabase + " SET "+AACompositioncolumn+"=%s WHERE UID =%s"
					mycursor.execute(UpdateStatement,InsertionValues)

					cnx.commit()
					
				except:
					pass
cnx.close()

print('Analysis Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-start_time))




