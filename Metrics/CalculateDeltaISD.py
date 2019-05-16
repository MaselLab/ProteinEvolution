import os, sys, mysql.connector, datetime

'''
Author : Sara Willis
Date   : May 13, 2019

This script is designed to calculate Delta ISD scores for entries in stored in MySQL Tables. 

Delta ISD is calculated as the difference in disorder prediction between a genuine protein and a version whose peptides have been randomly scrambled
'''




########################################################################
#                        User-Specific Information                     #                     
########################################################################
start_time = datetime.datetime.now()

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''


# Scrambled Tables -- where scrambled peptide sequences are stored with their output Tango scores
EnsemblScrambledTable = 'Genomes_Multicellular_ScrambledSequences' # Table Name
EnsemblScrambledISD = 'MeanISD_IUPred2_WithCys' 
EnsemblScrambledUID = 'ProteinTableUID' 


NCBIScrambledTable = 'NCBIGenomes_ScrambledSequences_Complete'
NCBIScrambledISD = 'MeanISD_IUPred2_WithCys'
NCBIScrambledUID = 'ProteinTableUID'

# Protein Metrics Tables -- The delta tango scores will be stored in these tables
EnsemblProteinMetricsTable = 'Genomes_Multicellular_ProteinMetrics'
EnsemblProteinMetricsISD = 'MeanISD_IUPred2_WithCys'
EnsemblProteinMetricsUID = 'ProteinTableUID'
EnsemblProteinMetricsDeltaISD = 'DeltaISD_IUPred2_WithCys'

NCBIProteinMetricsTable = 'NCBIGenomes_ProteinMetrics_Complete'
NCBIProteinMetricsISD = 'MeanISD_IUPred2_WithCys'
NCBIProteinMetricsUID = 'ProteinTableUID'
NCBIProteinMetricsDeltaISD = 'DeltaISD_IUPred2_WithCys'


# The user's input table information is stored as a dictionary for ease of use
DataTables = {'Ensembl':
              {'Scrambled':
               {'Table':EnsemblScrambledTable,
                'UID':EnsemblScrambledUID,
                'ISD': EnsemblScrambledISD},
                
               'ProteinMetrics':
               {'Table':EnsemblProteinMetricsTable,
                'UID':EnsemblProteinMetricsUID,
                'ISD':EnsemblProteinMetricsISD,
                'DeltaISD':EnsemblProteinMetricsDeltaISD}},
              
              'NCBI':
              {'Scrambled':
               {'Table': NCBIScrambledTable,
                'UID':NCBIScrambledUID,
                'ISD':NCBIScrambledISD},
               
               'ProteinMetrics':
               {'Table':NCBIProteinMetricsTable,
                'UID':NCBIProteinMetricsUID,
                'ISD': NCBIProteinMetricsISD,
                'DeltaISD':NCBIProteinMetricsDeltaISD}
              }}

########################################################################
#                        Program Executes Below                        #                     
########################################################################


# A connection to the MySQL database is established
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

for DataSource in DataTables:
    intermediate_time = datetime.datetime.now()

    # The data are extracted for each data source
    Scrambled = DataTables[DataSource]['Scrambled']['Table']
    ScrambledUID = DataTables[DataSource]['Scrambled']['UID']
    ScrambledISD = DataTables[DataSource]['Scrambled']['ISD']
               
    ProteinMetrics = DataTables[DataSource]['ProteinMetrics']['Table']
    ProteinMetricsUID = DataTables[DataSource]['ProteinMetrics']['UID']
    ProteinMetricsISD = DataTables[DataSource]['ProteinMetrics']['ISD']
    ProteinMetricsDeltaISD = DataTables[DataSource]['ProteinMetrics']['DeltaISD']
    
    # The relevant information is extracted 
    ExtractionStatement = "SELECT "+','.join([ProteinMetrics+"."+ProteinMetricsUID,ProteinMetrics+"."+ProteinMetricsISD,Scrambled+"."+ScrambledISD]) + " FROM "+Scrambled + " INNER JOIN " + ProteinMetrics + " ON " + ProteinMetrics+"."+ProteinMetricsUID +"="+Scrambled+"."+ScrambledUID
    
    print('Extracting Data From %s'%DataSource)
    sys.stdout.flush()
    
    mycursor.execute(ExtractionStatement)
    
    print('Data Extracted\nTime Taken: %s'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
    intermediate_time = datetime.datetime.now()

    # The extracted data are then parsed
    results = mycursor.fetchall()
    for result in results:
        UID,ISD_Real,ISD_Scrambled = result

        # Before any calculates are executed, the program makes sure values exist
        if ISD_Real != None and ISD_Scrambled != None:

            # Delta ISD is caluclated as the difference between the predicted ISD of a genuine protein
            # and the predicted ISD of it's scrambled version
            DeltaISD = float(ISD_Real)-float(ISD_Scrambled)

            # Once the delta ISD value has been calculated it's uploaded into the protein metrics table
            UpdateStatement = "UPDATE "+ProteinMetrics+" SET "+ProteinMetricsDeltaISD+"=%s "+" WHERE "+ProteinMetricsUID +"=%s"
            UpdateData = (DeltaISD,UID)
            mycursor.execute(UpdateStatement,UpdateData)
            cnx.commit()
    print('Delta ISD Values Updated\nTime Taken: %s'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
print('Script Complete\nTime Taken: %s'%(datetime.datetime.now()-start_time))
sys.stdout.flush()


cnx.close()
