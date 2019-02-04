import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib
from scipy import stats
import numpy as np
matplotlib.use('agg')
import matplotlib.pyplot as plt

'''
Author : Sara Willis
Date   : February 4, 2019

TO RUN:

        python3 <filename.py> --Metric [--y_min=VALUE --y_max=VALUE]

The purpose of this script is to generate scatterplots of the means and standard errors of a user-specified metric vs. Age for pfam domains across individual species. 

The user should specify the metric that should be plotted in the command line using the format specified under TO RUN. The options are the following:

   1) --ISD : The mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the domain

The user also has the option to make the various plots have uniform axes by specifying a y_min and y_max value for each figure. These are optional arguments. If none are specified, matplotlib uses the default settings. 
'''

# Enter MySQL connection information here
Database = ''
User = ''
Host = ''
Password = ''

# The database where the pfam UIDs, their ages, and their metrics are stored
PfamDataTable = 'PfamUIDsTable_EnsemblAndNCBI'

# These are the prefixes of the data tables where we will be selecting our data from. Data table name formatting: [Database]Genomes_DomainMetrics_Complete
SpeciesDatabases = ['Ensembl','NCBI']

# To customize the list of species you're generating your plots for, add a tuple to the dictionary below with the form (speciesUID, NewickSpeciesName).
# The tuple will need to go under the correct key pointing to which database the species was extracted from (NCBI or Ensembl) and must go inside the master list
def LoadSpeciesList(Database):
    SpeciesDictionary = {'NCBI':[(527,'Picoides_pubescens'),(392,'Oryza_sativa'),(366,'Capsicum_baccatum'),(464,'Aptenodytes_forsteri')],
                         'Ensembl':[(33,'Gallus_gallus'),(7,'BosTaurus'),(70,'Panthera_tigris'),(54,'Mus_musculus'),(14,'Cavia_porcellus'),(67,'Pan_paniscus'),(38,'Homo_sapiens'),(44,'Macaca_fascicularis'),(84,'Takifugu_rubripes'),(26,'Drosophila_melanogaster'),(209,'Drosophila_sechellia'),(283,'Aspergillus_clavatus'),(79,'Saccharomyces_cerevisiae')]}
    return SpeciesDictionary[Database]

# Metric Options
MetricOptions = {'ISD': ['MeanISD_IUPred2_WithCys', 'Mean ISD (IUPred 2)'],
                 'DensityOfAPRs': ['DensityOfAggProneRegions','Density of Aggregation-Prone Regions'],
                 'DensityOfAAsInAPRs': ['DensityOfAAsInAPRs',"Density of AAs in APRs"],
                 'Clustering_Trunc': ['NormalizedIndexOfDispersion_Trunc_FILVM', 'Clustering (Trunc,FILVM)'],
                 'Clustering_AllFrames': ['NormalizedIndexOfDispersion_AllPhases_FILVM','Clustering (All Frames, FILVM)'],
                 'Length': ['DomainLength', 'Domain Length']}


##########################################################################

# The program begins by trying to determine which metric the user would like plotted vs. time
# If the metric is located, the program proceeds. Otherwise, the program exits and notifies the
# user that there was an issue.

# The program also allows the user to input max and min y-values so the plots across species share axes to make interpretation easier

# The program checks to see if there are any user-arguments. If there aren't, the program hasn't been told which metric to plot. In this
# case, the program notifies the user that it cannot run and shuts down.
Input = sys.argv
if len(Input) == 1:
    print('\n\n\nInvalid format\n\nTo generated figures, a metric must be selected\n\nTo run, use the following format:\n\n\t' + '\33[34m' + 'python3 script_name.py --metric [--y_min=VALUE --y_max=VALUE]'+'\33[0m'+'\n\nFor help, use --options\n\n\n')
    sys.exit(0)
    
# The user can request help from the command line using the option 'options'. This will provide a list of the metrics the user can
# specify for the program to plot
if '--options' in Input:
    print('\n\nThe following metrics options are available to plot vs. age:\n\n\t--ISD\n\t--DensityOfAPRs\n\t--DensityOfAAsInAPRs\n\t--Clustering_Trunc\n\t--Clustering_AllFrames\n\t--Length\n\nIf the user wishes to specify min/max values for the y-axis to force\nthe species plots to be uniformly scaled, use the arguments:\n\n\t--y_min=VALUE\n\t--y_max=VALUE\n\n\n')
    sys.exit(0)

MetricOption = None
y_min = None
y_max = None

# If there are user options specified, the program checks to see if the metric has been specified and whether there are y_min and y_max values specified 
for inputItem in Input:
    if inputItem.replace('--','') in MetricOptions:
        MetricOption = inputItem.replace('--','')
    if 'y_min' in inputItem:
        try:
            y_min = float(inputItem.split('=')[1])
        except:
            # In the case that the user has input non-numeric values for either y_min or y_max, the program ignores
            # the user-specifications and notifies the user that it cannot use the input value(s)
            print('\n\n'+'\33[31m'+'Warning: Not a valid y_min value\nValue must be numeric. Ignoring.'+'\33[0m'+'\n\n')
            pass
    if 'y_max' in inputItem:
        try:
            y_max = float(inputItem.split('=')[1])
        except:
            print('\n\n'+'\33[31m'+'Warning: Not a valid y_max value.\nValue must be numeric. Ignoring.'+'\33[0m'+'\n\n')
            pass
try:
    # If the metric was not in the dictionary MetricOptions, something has gone wrong and MetricOption remains None.
    # If this is the case, the user is notified that the metric option is not valid and the program exits
    MetricsColumn = MetricOptions[MetricOption.replace('-','')][0]
    MetricsLabel = MetricOptions[MetricOption.replace('-','')][1]
except:
    print('\n\nInvalid metric option selected\n\nPlease specify a valid metric\n\nFor help, use --options\n\n')
    sys.exit(0)

# A connection to the MySQL database is established
startTime = datetime.datetime.now()
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We'll search both the Ensembl and the NCBI data tables 
for SpeciesDatabase in SpeciesDatabases:
    print('\nGenerating Figures!\n\n')
    speciesList = LoadSpeciesList(SpeciesDatabase)
    # Each species specified by the user will get its own unique plot
    for entry in speciesList:
        SpeciesUID = entry[0]
        SpeciesName = entry[1]
        # The species UID is used to extract species-specific data points. Species UID should be
        # listed as an index in the data table, otherwise this process could take a very long time!
        ExtractionStatement = "SELECT PfamUID,PassedIUPredFilters," + MetricsColumn + " FROM "+SpeciesDatabase+"Genomes_DomainMetrics_Complete WHERE SpeciesUID = %s"%SpeciesUID
        print('\n\nExtracting data for %s'%SpeciesName)
        mycursor.execute(ExtractionStatement)
        results = mycursor.fetchall()

        # Because Pfams may show up multiple times in a species, we will create a dictionary with
        # an entry for each Pfam UID where we will store the metric value associated with that
        # particular pfam instance. Once the Pfam dictionary is compiled, the average value of all
        # metrics associated with that Pfam UID is calculated. 
        PfamDictionary = {}
        
        for result in results:
            if result != None:
                PfamUID,PassedFilters,Metric = result
                # We only want to collect metrics from sequences that have passed a preliminary quality control
                # PassedFilters indicates whether a sequence was of sufficent quality to be considered. To pass
                # the quality filter, the coding sequence needs to meet all three of the following requirements:
                #   1) Starts with a start codon
                #   2) Contains no in-frame stop codons
                #   3) Is a multiple of three in length
                if Metric != None and PassedFilters != int(False):
                    Metric = float(Metric)
                    if PfamUID not in PfamDictionary:
                        PfamDictionary[PfamUID] = [Metric]
                    else:
                        PfamDictionary[PfamUID].append(Metric)
        print('Extraction Complete\nTime Taken: %s\n\nCreating Dictionary\n\n'%(datetime.datetime.now()-startTime))
        currentTime = datetime.datetime.now()
        # Once the the pfam database is compiled, the pfams that share common ages are consolidated into a dictionary.
        # Once the dictionary is complete, the average of the pfam values sharing that age is found. This is the value
        # that is ultimately plotted
        AgeDictionary = {}

        # To assign ages to the pfams, we need to access the database where that information is stored
        for pfam in PfamDictionary:
            SelectPfamAgeStatement = "SELECT Age_MY FROM "+PfamDataTable+" WHERE PfamUID = '%s'"%pfam
            mycursor.execute(SelectPfamAgeStatement)
            age = mycursor.fetchone()[0]
            if age not in AgeDictionary:
                AgeDictionary[age] = [np.mean(PfamDictionary[pfam])]
            else:
                AgeDictionary[age].append(np.mean(PfamDictionary[pfam]))

        # Once the age dictionary is complete, the program plots the figures for us
        # It assembles a variety of lists to plot the data
        print('Data Extracted\nTime Taken: %s\n\nCreating Figure'%(datetime.datetime.now()-currentTime))

        # The first set of lists store mean metric values (y), age (x), and standard error (e)
        x = []
        y = []
        e = []

        # If there's only one metric value associated with a particular age, then no standard error value
        # is available, so we store these values in a separate set of lists so they can be plotted using a
        # slightly different format
        x_1 = []
        y_1 = []

        # Constructs lists
        for age in AgeDictionary:
            if len(AgeDictionary[age]) != 1:
                standardError = stats.sem(AgeDictionary[age])
                x.append(age)
                y.append(np.mean(AgeDictionary[age]))
                e.append(standardError)
            else:
                x_1.append(age)
                y_1.append(AgeDictionary[age][0])

        # The scatterplot + errorbar plotting function plt.errorbar can be used
        # for entries with an associated error term. It is not permitted to not have
        # error terms.
        plt.errorbar(x, y, e, linestyle='None', markersize=5,marker='o')
        # Points with no associated error term are plotted using the plt.scatter function
        # and are plotted in a slightly paler color using a smaller marker
        plt.scatter(x_1, y_1, c = '#3FC0DA' ,linestyle = 'None', marker='.')
        plt.xlabel('Age (MY)')
        plt.ylabel(MetricsLabel)
        if y_max != None and y_min == None:
            plt.ylim(0,y_max)
        if y_max != None and y_min != None:
            plt.ylim(y_min,y_max)
        plt.xlim(0,2300)
        plt.title('Species: %s'%SpeciesName)
        # Because plots are generated using this script for multiple species, and different runs
        # may produce many different metrics plots, generating all figures to the same directory
        # may be untenable. To try to reduce the madness, the program writes directories to store
        # results based on the metrics that are being plotted. 
        if os.path.isdir("./%s_PfamOnly"%MetricOption) == False:
            os.system('mkdir %s_PfamOnly'%MetricOption)
        plt.savefig("./%s_PfamOnly/%s"%(MetricOption,SpeciesName))
        plt.close()
              
print("Done!\nTime Taken: %s"%(datetime.datetime.now()-startTime))

        
