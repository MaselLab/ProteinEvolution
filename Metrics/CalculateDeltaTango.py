import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib, requests

'''
Author : Sara Willis
Date   : May 13, 2019

This script is designed to calculate Delta Tango values for entries in stored in MySQL Tables. 

Delta tango calculates two scores:
    
    1) Density of APRs (original protein) - Density of APRs (scrambled protein)

    2) Density of AAs in APRs (original protein) - Density of APRs (scrambled protein)

The delta tango values are then stored the defined protein metrics table

'''




########################################################################
#                        User-Specific Information                     #                     
########################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''


# Scrambled Tables -- where scrambled peptide sequences are stored with their output Tango scores
EnsemblScrambledTable = 'Genomes_Multicellular_ScrambledSequences' # Table Name
EnsemblScrambledDAPRs = 'DensityOfAggProneRegions' 
EnsemblScrambledDAAAPRs = 'DensityOfAAsInAPRs'
EnsemblScrambledUID = 'ProteinTableUID' 


NCBIScrambledTable = 'NCBIGenomes_ScrambledSequences_Complete'
NCBIScrambledDAPRs = 'DensityOfAggProneRegions'
NCBIScrambledDAAAPRs = 'DensityOfAAsInAPRs'
NCBIScrambledUID = 'ProteinTableUID'

# Protein Metrics Tables -- The delta tango scores will be stored in these tables
EnsemblProteinMetricsTable = 'Genomes_Multicellular_ProteinMetrics'
EnsemblProteinMetricsDAPRs = 'DensityOfAggProneRegions'
EnsemblProteinMetricsDAAAPRs = 'DensityOfAAsInAPRs'
EnsemblProteinMetricsUID = 'ProteinTableUID'
EnsemblProteinMetricsDeltaDAPRs = 'DeltaTango_DensityOfAggProneRegions'
EnsemblProteinMetricsDeltaDAAAPRs = 'DeltaTango_DensityOfAAsInAPRs'

NCBIProteinMetricsTable = 'NCBIGenomes_ProteinMetrics_Complete'
NCBIProteinMetricsDAPRs = 'DensityOfAggProneRegions'
NCBIProteinMetricsDAAAPRs = 'DensityOfAAsInAPRs'
NCBIProteinMetricsUID = 'ProteinTableUID'
NCBIProteinMetricsDeltaDAPRs = 'DeltaTango_DensityOfAggProneRegions'
NCBIProteinMetricsDeltaDAAAPRs = 'DeltaTango_DensityOfAAsInAPRs'

# The user's data are stored in a dictionary so they're easier to parse
DataTables = {'Ensembl':
              {'Scrambled':
               {'Table':EnsemblScrambledTable,
                'UID':EnsemblScrambledUID,
                'DAPRs': EnsemblScrambledDAPRs,
                'DAAAPRs':EnsemblScrambledDAAAPRs},
                
               'ProteinMetrics':
               {'Table':EnsemblProteinMetricsTable,
                'UID':EnsemblProteinMetricsUID,
                'DAPRs':EnsemblProteinMetricsDAPRs,
                'DAAAPRs':EnsemblProteinMetricsDAAAPRs,
                'DeltaDAPRs':EnsemblProteinMetricsDeltaDAPRs,
                'DeltaDAAAPRs':EnsemblProteinMetricsDeltaDAAAPRs}},
              
              'NCBI':
              {'Scrambled':
               {'Table': NCBIScrambledTable,
                'UID':NCBIScrambledUID,
                'DAPRs':NCBIScrambledDAPRs,
                'DAAAPRs':NCBIScrambledDAAAPRs},
               
               'ProteinMetrics':
               {'Table':NCBIProteinMetricsTable,
                'UID':NCBIProteinMetricsUID,
                'DAPRs': NCBIProteinMetricsDAPRs,
                'DAAAPRs':NCBIProteinMetricsDAAAPRs,
                'DeltaDAPRs':NCBIProteinMetricsDeltaDAPRs,
                'DeltaDAAAPRs':EnsemblProteinMetricsDeltaDAAAPRs}
              }}

########################################################################
#                        Program Executes Below                        #                     
########################################################################

start_time = datetime.datetime.now()
print('Program Executing\nCurrent Time: %s'%start_time)
sys.stdout.flush()

# A connection is established to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Each data source's delta tango values will be calculated (either Ensembl or NCBI)
for DataSource in DataTables:
    intermediate_time = datetime.datetime.now()
    # The table and column names are loaded
    Scrambled = DataTables[DataSource]['Scrambled']['Table']
    ScrambledUID = DataTables[DataSource]['Scrambled']['UID']
    ScrambledDAPRs = DataTables[DataSource]['Scrambled']['DAPRs']
    ScrambledDAAAPRs = DataTables[DataSource]['Scrambled']['DAAAPRs']
               
    ProteinMetrics = DataTables[DataSource]['ProteinMetrics']['Table']
    ProteinMetricsUID = DataTables[DataSource]['ProteinMetrics']['UID']
    ProteinMetricsDAPRs = DataTables[DataSource]['ProteinMetrics']['DAPRs']
    ProteinMetricsDAAAPRs = DataTables[DataSource]['ProteinMetrics']['DAAAPRs']

    ProteinMetricsDeltaDAPRs = DataTables[DataSource]['ProteinMetrics']['DeltaDAPRs']
    ProteinMetricsDeltaDAAAPRs = DataTables[DataSource]['ProteinMetrics']['DeltaDAAAPRs']

    # and an extraction statement is defined and executed, extracting the relevant information
    # from MySQL
    ExtractionStatement = "SELECT "+','.join([ProteinMetrics+"."+ProteinMetricsUID,ProteinMetrics+"."+ProteinMetricsDAPRs,ProteinMetrics+"."+ProteinMetricsDAAAPRs,Scrambled+"."+ScrambledDAPRs,Scrambled+"."+ScrambledDAAAPRs]) + " FROM "+Scrambled + " INNER JOIN " + ProteinMetrics + " ON " + ProteinMetrics+"."+ProteinMetricsUID +"="+Scrambled+"."+ScrambledUID
    print('Extracting Data From %s\n'%(datetime.datetime.now()-start_time))
    sys.stdout.flush()
    mycursor.execute(ExtractionStatement)

    print('Data Extracted\nTime Taken: %s\n\nCalculating Delta Tango Values\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
    intermediate_time = datetime.datetime.now()
    # The results are then collected and sorted through
    results = mycursor.fetchall()
    for result in results:
        UID,DAPRs_Real,DAAAPRs_Real,DAPRs_Scrambled,DAAAPRs_Scrambled = result
        TangoValues = [DAPRs_Real,DAPRs_Scrambled,DAAAPRs_Real,DAAAPRs_Scrambled]
        
        # Before attempting to perform any operations on the Tango values, we make sure
        # they all exist
        if None not in TangoValues:
            DeltaDAPRs = float(DAPRs_Real)-float(DAPRs_Scrambled)
            DeltaDAAAPRs = float(DAAAPRs_Real)-float(DAPRs_Scrambled)

            # The delta tango values are then inserted back into the protein metrics table
            UpdateStatement = "UPDATE "+ProteinMetrics+" SET "+ProteinMetricsDeltaDAPRs+"=%s"+","+ProteinMetricsDeltaDAAAPRs+"=%s"+" WHERE "+ProteinMetricsUID +"=%s"
            UpdateData = (DeltaDAPRs,DeltaDAAAPRs,UID)
            mycursor.execute(UpdateStatement,UpdateData)
            cnx.commit()
    print('Delta Tango Values Uploaded\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
print('\nProgram Complete\nTime Taken: %s'%(datetime.datetime.now()-start_time))
sys.stdout.flush()



cnx.close()
