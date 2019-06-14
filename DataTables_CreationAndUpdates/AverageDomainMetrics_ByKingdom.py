import os, sys, mysql.connector, datetime
import numpy as np

'''
Author : Sara Willis
Date   : June 13, 2019

This script is designed to calculate the means and variances of a selected metric for domains by kingdom. It does this by extracting all instances of each pfam identified as belonging to a particular kingdom. The mean and variance of the selected metric of that particular domain is then taken over all occurences within that kingdom and the results are uploaded to a MySQL table. 

The current available metrics are ISD and Clustering. They can either be run individually or together.

** Note -- this code can be optimized so that both metrics are extracted using a single MySQL query rather than using a nest for loop. It was written this way as the option to include multiple metrics was included as an afterthought. It hasn't been optimized as of yet because the full script runs at ~1.5 minutes per metric and so is not of particularly high priority relative to other tasks. Still, it bothers me and so I'm including this comment.

'''
##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################


# Either 'Clustering', 'ISD', or both as comma-delimited strings 
Metrics = 'Clustering','ISD'

# User-Specific SQL Connection Information
Database                             = ''
User                                 = ''
Host                                 = ''
Password                             = ''

# Domain Metrics Table
DomainMetricTables                   = 'NCBIGenomes_DomainMetrics_Complete','Genomes_Multicellular_DomainMetrics'
SpeciesUID_DomainTable               = 'SpeciesUID'
PfamUID_DomainTable                  = 'PfamUID'
ISD_DomainTable                      = 'MeanISD_IUPred2_WithCys'
Clustering_DomainTable               = 'NormalizedIndexOfDispersion_Trunc_FILVM'
DomainMetricTables_ProteinUID        = 'ProteinTableUID'

# Unique Protein Transcripts Table
UniqueTranscriptsTables              = 'UniqueProteinGenePairs_Table_NCBI','UniqueProteinGenePairs_Table_Ensembl'
UniqueTranscriptUIDs                 = 'UID'

# Pfam UIDs Table -- Where all processed data are stored
PfamUIDsTable                        = 'PfamUIDsTable_EnsemblAndNCBI'
PfamUID_PfamTable                    = 'PfamUID'
ISD_Animal_PfamTable                 = 'MeanISD_AnimalSpecific'
ISD_Variance_Animal_PfamTable        = 'ISDVariance_AnimalSpecific'
ISD_Plant_PfamTable                  = 'MeanISD_PlantSpecific'
ISD_Variance_Plant_PfamTable         = 'ISDVariance_PlantSpecific'

Clustering_Animal_PfamTable          = 'MeanClustering_AnimalSpecific'
Clustering_Variance_Animal_PfamTable = 'ClusteringVariance_AnimalSpecific'
Clustering_Plant_PfamTable           = 'MeanClustering_PlantSpecific'
Clustering_Variance_Plant_PfamTable  = 'ClusteringVariance_PlantSpecific'

# Species Table -- Used to determine pfam species' kingdom
SpeciesTable                         = 'SpeciesList'
SpeciesUID_SpeciesTable              = 'SpeciesUID'
Kingdom_SpeciesTable                 = 'Kingdom'




##########################################################################################
#                                Program Executes Below                                  #
##########################################################################################



start_time = datetime.datetime.now()

# Establish a connection to the mysql database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)


# If the user has specified more than one domain table, the source tables variable will be
# a tuple, so the domains tables and their corresponding unique transcript tables are
# zipped together. If only one table has been specified, it is converted into a list object
# with its unique transcripts table so it can be manipulated in the same way as when two or
# more tables are specified
if type(DomainMetricTables) == tuple:
    DataTables = list(zip(DomainMetricTables,UniqueTranscriptsTables))
else:
    DataTables = [DomainMetricTables,UniqueTranscriptsTables]

# More than one metric can be supplied by the user. The metrics variable is dealt with in the
# same way as the DataTables. This could really be optimized so that there are different
# MySQL extraction statements so that everything is pulled at once rather than using nested
# for loops. This, however, is not a high priority since this script takes ~1.5 minutes per
# data table. Other tasks beckon!
if type(Metrics) != tuple:
    Metrics = [Metrics]

# The script makes sure that the metrics that have been supplied by the user are usable
if set(Metrics).issubset({'ISD','Clustering'}) == False:
    print(5*'\n'+40*'_'+"\nMetric Option Not Recognized\n\nMetric Options:\n\n\t* 'ISD'\n\t* 'Clustering'\n"+40*'_'+5*'\n')
    sys.stdout.flush()
    cnx.close()
    sys.exit(0)

# Here is where optimization could occur (were someone so inclined...). This for loop could
# be removed and a customizable SQL extraction statement could be made
for Metric in Metrics:
    print('\nExtracting %s Values\n'%Metric)
    sys.stdout.flush()

    # Each instance of a Pfam is stored as an entry in a dictionary with its metric values
    # stored in a list. Following the sorting process, the mean and variance of that metric
    # for that particular pfam is found
    PfamUIDs = {}

    # We label the variables appropriately so the correct information is extracted from the
    # MySQL table
    if Metric == 'Clustering':
        MetricColumn = Clustering_DomainTable
        Plant_Mean = Clustering_Plant_PfamTable
        Plant_Variance = Clustering_Variance_Plant_PfamTable
        Animal_Mean = Clustering_Animal_PfamTable
        Animal_Variance = Clustering_Variance_Animal_PfamTable

    elif Metric == 'ISD':
        MetricColumn = ISD_DomainTable
        Plant_Mean =ISD_Plant_PfamTable
        Plant_Variance =ISD_Variance_Plant_PfamTable
        Animal_Mean =ISD_Animal_PfamTable
        Animal_Variance = ISD_Variance_Animal_PfamTable

    # Because the data are split up between multiple tables, we need to configure an
    # an extraction statement for each source data table and will need to pull all entries
    # and store the output before we calculate means and variances
    for Table in DataTables:
        DomainMetricTable, UniqueTable = Table
        intermediate_time = datetime.datetime.now()
        print('Extracting Data From: %s\n'%DomainMetricTable)
        sys.stdout.flush()

        # Data extracted: PfamUID, Metric and Kingdom -- these are pulled together by linking
        # the species UID in the domain metrics table with the kingdom identified in the
        # species table. Only one transcript per protein is chosen using the "unique" table
        PfamUIDsExtractionStatement = "SELECT "+','.join([DomainMetricTable+'.'+i for i in [PfamUID_DomainTable,MetricColumn]]+[SpeciesTable+'.'+i for i in [Kingdom_SpeciesTable]])+" FROM "+DomainMetricTable+" INNER JOIN "+SpeciesTable+" ON %s.%s"%(DomainMetricTable,SpeciesUID_DomainTable)+"=%s.%s"%(SpeciesTable,SpeciesUID_SpeciesTable) + " INNER JOIN "+UniqueTable+" ON %s.%s"%(DomainMetricTable,DomainMetricTables_ProteinUID)+"=%s.%s"%(UniqueTable,UniqueTranscriptUIDs) + " WHERE "+DomainMetricTable+".UID<100"
        mycursor.execute(PfamUIDsExtractionStatement)
        
        print('Data Successfully Extracted\nTime Taken: %s\n\nSorting Entries'%(datetime.datetime.now()-intermediate_time))
        sys.stdout.flush()
        
        intermediate_time = datetime.datetime.now()
        
        results = mycursor.fetchall()
        
        for result in results:
            PfamUID, ExtractedMetric, Kingdom = result
            
            # Sometimes things go wrong and a particular metric value is not available.
            # We ignore these cases
            if ExtractedMetric != None:

                # The kingdoms we're interested in are plants and animals so we clump them
                # into individual categories in our dictionary
                if Kingdom == 'plant':
                    if PfamUID not in PfamUIDs:
                        PfamUIDs[PfamUID] = {'Plant': [ExtractedMetric],
                                             'Animal':[]}
                    else:
                        PfamUIDs[PfamUID]['Plant'].append(ExtractedMetric)

                # Animals are not stored with the identifier 'Animal' in the database but are
                # partitioned into either 'vertebrate' or 'invertebrate'. We concateneate these
                # to get the desired output
                elif Kingdom == 'invertebrate' or Kingdom == 'vertebrate':
                    if PfamUID not in PfamUIDs:
                        PfamUIDs[PfamUID] = {'Plant': [],
                                             'Animal':[ExtractedMetric]}
                    else:
                        PfamUIDs[PfamUID]['Animal'].append(ExtractedMetric)
                        
        print('Results Sorted\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
        sys.stdout.flush()
        
    intermediate_time = datetime.datetime.now()
    
    print('All Entries Sorted\nTaking Averages and Uploading Results\n')
    sys.stdout.flush()
    
    for PfamUID in PfamUIDs:
        # Not all Pfams are present in particular kingdoms so when we're sorting through all Pfams
        # we need to keep that in mind and not try to find the average value of an empty list
        PlantMean = np.mean(PfamUIDs[PfamUID]['Plant']) if len(PfamUIDs[PfamUID]['Plant']) != 0 else None
        PlantVariance = np.var(PfamUIDs[PfamUID]['Plant']) if len(PfamUIDs[PfamUID]['Plant']) != 0 else None
        AnimalMean = np.mean(PfamUIDs[PfamUID]['Animal']) if len(PfamUIDs[PfamUID]['Animal'])!= 0 else None
        AnimalVariance = np.var(PfamUIDs[PfamUID]['Animal']) if len(PfamUIDs[PfamUID]['Animal']) != 0 else None
        
        # The results are then patched up with a SQL extraction statement and uploaded to the Pfam UIDs table
        UpdateStatement = "UPDATE "+PfamUIDsTable+" SET "+Animal_Mean+"=%s,"+Plant_Mean+"=%s,"+Animal_Variance+"=%s,"+Plant_Variance+"=%s WHERE "+PfamUID_PfamTable+"=%s"
        UpdateTuple = (AnimalMean,PlantMean,AnimalVariance,PlantVariance,PfamUID)
        mycursor.execute(UpdateStatement,UpdateTuple)
        cnx.commit()
    
print('Averages Taken and Uploaded\nTime Taken: %s\n\nScript Complete!\nTotal Time Taken: %s'%(datetime.datetime.now()-intermediate_time,datetime.datetime.now()-start_time))
sys.stdout.flush()
                
cnx.close()
