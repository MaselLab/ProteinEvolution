import os, json, sys, csv, mysql.connector, datetime, copy
from scipy import stats
import numpy as np


'''
Author : Sara Willis
Date   : February 19, 2019


This script is written to create a datatable that has a unique entry for each Pfam UID that appears in the Ensembl and NCBI datatables. 

It accomplishes this by extracting every domain entry in the datatables EnsemblGenomes_DomainMetrics_Complete and NCBIGenomes_DomainMetrics_Complete. In assembles these entries into a dictionary so that the metrics for each occurence of each PfamUID are stored in a list associated with that Pfam domain. 


'''

###########################################################################################
#                                User-Specific Information                                #        
###########################################################################################

# User's MySQL Connection Information 
Database = ''
User = ''
Host = ''
Password = ''

EnsemblDomainDataTable = 'EnsemblGenomes_DomainMetrics_Complete'
NCBIDomainDataTable = 'NCBIGenomes_DomainMetrics_Complete'

SpeciesDataTable = 'SpeciesList'

Verbose = True

###########################################################################################
#                             Program Executes Below                                      #         
###########################################################################################

startTime = datetime.datetime.now()
# Establishes a connection to the MySQL database 
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

DataTables = [EnsemblDomainDataTable,NCBIDomainDatatTable]

PfamDictionary = {}
SpeciesDictionary = {}

# Searches both data tables for all pfams that show up
for DataTable in DataTables:
    currentTime = datetime.datetime.now()
    if Verbose == True:
        print('\n\nExtracting from %s\n\n\n'%DataTable)
        sys.stdout.flush()
        
    # Extracts entries one at a time
    mycursor.execute("SELECT MAX(UID) FROM "+DataTable)
    MaxUID = mycursor.fetchone()[0]
    currentTime = datetime.datetime.now()
    
    for uid in range(1,MaxUID):
        
        if uid%1000000 == 0 and Verbose == True:
            print('%s processed'%uid,end = '\r')
            sys.stdout.flush()
        # Selects PfamUIDs along with all their metrics
        ExtractionStatement = "SELECT PfamUID,SpeciesUID,DomainLength,MeanISD_IUPred2_WithCys,DensityOfAggProneRegions,DensityOfAAsInAPRs,NormalizedIndexOfDispersion_Trunc_FILVM,NormalizedIndexOfDispersion_AllPhases_FILVM,PassedIUPredFilters FROM "+DataTable+" WHERE UID = %s"%uid
        mycursor.execute(ExtractionStatement)

        mycursor.execute(ExtractionStatement)
        result = mycursor.fetchone()
        if result != None:
            PfamUID,SpeciesUID,Length,ISD,DensityOfAPRs,DensityOfAAsInAPRs,Clustering_Trunc,Clustering_AllPhases,PassedIUPredFilters = result

            # We then go through the process of creating a pfam dictionary. Each key is a unique Pfam ID and points to a nested dictionary
            # that contains its metrics and the species that it shows up in 
            if PfamUID not in PfamDictionary:
                if PassedIUPredFilters == int(True):
                    PfamDictionary[PfamUID] = {'SpeciesUIDs':set([SpeciesUID]),
                                               'ISD':[ISD],
                                               'Length':[Length],
                                               'DensityOfAPRs':[DensityOfAPRs],
                                               'DensityOfAAsInAPRs':[DensityOfAAsInAPRs],
                                               'Clustering_Trunc':[Clustering_Trunc],
                                               'Clustering_AllPhases':[Clustering_AllPhases]}
                else:
                    PfamDictionary[PfamUID] = {'SpeciesUIDs':set([SpeciesUID]),
                                               'ISD':[],
                                               'Length':[],
                                               'DensityOfAPRs':[],
                                               'DensityOfAAsInAPRs':[],
                                               'Clustering_Trunc':[],
                                               'Clustering_AllPhases':[]}

            else:
                if PassedIUPredFilters == int(True):
                    PfamDictionary[PfamUID]['SpeciesUIDs'].add(SpeciesUID)
                    PfamDictionary[PfamUID]['ISD'].append(ISD)
                    PfamDictionary[PfamUID]['Length'].append(Length)
                    PfamDictionary[PfamUID]['DensityOfAPRs'].append(DensityOfAPRs)
                    PfamDictionary[PfamUID]['DensityOfAAsInAPRs'].append(DensityOfAAsInAPRs)
                    PfamDictionary[PfamUID]['Clustering_Trunc'].append(Clustering_Trunc)
                    PfamDictionary[PfamUID]['Clustering_AllPhases'].append(Clustering_AllPhases)
                else:
                    PfamDictionary[PfamUID]['SpeciesUIDs'].add(SpeciesUID)
    if Verbose == True:
        print('Data Extracted From %s\nTotal Time Elapsed: %s\n\n'%(DataTable,datetime.datetime.now()-startTime))
        sys.stdout.flush()

# We want the names of the species the Pfams show up in, not the species UIDs, so we go through the species database and
# make a dictionary where each species UID points to the species' newick name
mycursor.execute('SELECT SpeciesUID,NewickSpeciesName FROM '+SpeciesDataTable)
results = mycursor.fetchall()
for result in results:
    SpeciesUID, NewickSpecies = result
    SpeciesDictionary[SpeciesUID] = NewickSpecies


# We then add the names of the species to each pfam as well as finding the mean metrics associated with that pfam
for PfamUID in PfamDictionary:
    speciesList = []
    for speciesUID in PfamDictionary[PfamUID]['SpeciesUIDs']:
        speciesList.append(SpeciesDictionary[speciesUID])
    speciesString = ','.join(speciesList)
    MeanISD = np.mean([i for i in PfamDictionary[PfamUID]['ISD'] if i != None])
    try:
        ISDVariance = np.var([i for i in PfamDictionary[PfamUID]['ISD'] if i != None])
    except:
        ISDVariance = None
    MeanLength = np.mean([i for i in PfamDictionary[PfamUID]['Length'] if i != None])
    MeanDensityOfAPRs = np.mean([i for i in PfamDictionary[PfamUID]['DensityOfAPRs'] if i != None])
    MeanDensityOfAAsInAPRs = np.mean([i for i in PfamDictionary[PfamUID]['DensityOfAAsInAPRs'] if i != None])
    MeanClusteringTrunc = np.mean([i for i in PfamDictionary[PfamUID]['Clustering_Trunc'] if i != None])
    MeanClusteringAllPhases = np.mean([i for i in PfamDictionary[PfamUID]['Clustering_AllPhases'] if i != None])
    NumberOfAssociatedSpecies = len(speciesList)
    InsertionStatement = "INSERT INTO PfamUIDsTable_EnsemblAndNCBI_copy_copy (PfamUID,Species,NumberOfAssociatedSpecies,ISDOverAllSpecies_Mean,ISDOverAllSpecies_Variance,ClusteringOverAllSpecies_Trunc_Mean,ClusteringOverAllSpecies_AllFrames_Mean,Aggregation_Tango_DensityOfAPRs_Mean,Aggregation_Tango_DensityOfAAsInAPRs_Mean,MeanDomainLength) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    InsertionData = (PfamUID, speciesString,NumberOfAssociatedSpecies,MeanISD,ISDVariance,MeanClusteringTrunc,MeanClusteringAllPhases,MeanDensityOfAPRs,MeanDensityOfAAsInAPRs,MeanLength)
    ReformattedInsertionData = tuple([None if str(i) == 'nan' else i for i in InsertionData])
    try:
        mycursor.execute(InsertionStatement,ReformattedInsertionData)
        cnx.commit()
    except Exception as e:
        if Verbose == True:
            print(e)
        



cnx.close()
