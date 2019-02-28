import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib
from scipy import stats
import numpy as np
matplotlib.use('agg')
import matplotlib.pyplot as plt


'''
Author : Sara Willis
Date   : February 1, 2019
-------------------------

* Note -- there is homology group work that's being done that will require this script to be updated


TO RUN: 
        python3 <filename.py> --Metric

Choose metrics option from below or type --options


The purpose of this script is to plot the means and standard errors of various protein metrics vs. the age of that protein.

The metrics determined for each protein are for that entire protein and not for any specific domains contained within it. 

There are a variety of metrics that can be plotted. The desired metric is specified by the user as a command-line argument. They are the following:

   1) --ISD : The mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the amino acid

'''


##########################################################################
#                               User Data                                #
##########################################################################


# SQL Connection Information 
Database = ''
User = ''
Host = ''
Password = ''


# If the user needs to make modifications to the names of the columns where the metrics are stored,
# or wishes to change the labeling scheme of their plot, they can make modifications to the dictionary
# below. The general format is:
#    Key = Input option
#    List entry 1: SQL column Name
#    List entry 2: Metric label on plot
MetricOptions = {'ISD': ['MeanISD_IUPred2_WithCys', 'Mean ISD (IUPred 2)'],
                 'DensityOfAPRs': ['DensityOfAggProneRegions','Density of Aggregation-Prone Regions'],
                 'DensityOfAAsInAPRs': ['DensityOfAAsInAPRs',"Density of AAs in APRs"],
                 'Clustering_Trunc': ['NormalizedIndexOfDispersion_Trunc_FILVM', 'Clustering (Trunc,FILVM)'],
                 'Clustering_AllFrames': ['NormalizedIndexOfDispersion_AllPhases_FILVM','Clustering (All Frames, FILVM)'],
                 'Length': ['ProteinLength', 'Protein Length']}




##########################################################################
#                           Program Executes Below                       #
##########################################################################


# The program begins by trying to determine which metric the user would like plotted vs. time
# If the metric is located, the program proceeds. Otherwise, the program exits and notifies the
# user that there was an issue
try:
    Metric = sys.argv[1]
except:
    print('\n------------------\n\nInvalid format\n\nTo generated figures, a metric must be selected\n\nTo run, use the following format:\n\n' + '\x1b[6;30;43m' + '<python3 command> <script name>.py --[metric option]'+'\x1b[0m'+'\n\nFor help, use --options\n\n------------------\n')
    sys.exit(0)

if Metric == '--options':
    print('\n\nThe following options are available to plot vs. age:\n\n--ISD\n--DensityOfAPRs\n--DensityOfAAsInAPRs\n--Clustering_Trunc\n--Clustering_AllFrames\n--Length\n\n')
    sys.exit(0)
try:
    MetricsColumn = MetricOptions[Metric.replace('-','')][0]
    MetricsLabel = MetricOptions[Metric.replace('-','')][1]
except:
    #print('\n\nInvalid metric option selected\n\nPlease use correct format:\n'+'\x1b[6;30;42m' +'<python3 command> ScriptName.py --option\n\n'+'\x1b[0m'+ '\nPlease specify a metric\n\nFor help, use --options\n\n')
    print('\n\nInvalid metric option selected\n\nPlease specify a valid metric\n\nFor help, use --options\n\n')
    sys.exit(0)

startTime = datetime.datetime.now()

# The script then establishes a connection with the SQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Each occurence of a metric is stored with all others of its age. The average value is found once all
# data entries are stored
AgeDictionary = {}

# The program extracts all data values from the Ensembl and NCBI databases
SpeciesDatabases = ['Ensembl','NCBI']
for SpeciesDatabase in SpeciesDatabases:
    # We want to not only extract the ages and the metric, but also a quality control.
    # The quality control PassedIUPredFilters is an indication that a coding sequence
    # was low quality. If this is the case, we ignore the protein
    ExtractionStatement = "SELECT AgeOfOldestPfam,PassedIUPredFilters," + MetricsColumn + " FROM "+SpeciesDatabase+"Genomes_ProteinMetrics_Complete "
    print('\n\nExtracting data from %s\n\n'%SpeciesDatabase)
    mycursor.execute(ExtractionStatement)
    results = mycursor.fetchall()
    for result in results:
        if result != None:
            Age,PassedFilters,ExtractedMetric = result
            if ExtractedMetric != None and PassedFilters != int(False):
                ExtractedMetric = float(ExtractedMetric)
                if Age not in AgeDictionary:
                    AgeDictionary[Age] = [ExtractedMetric]
                else:
                    AgeDictionary[Age].append(ExtractedMetric)
    print('Extraction Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-startTime))
    currentTime = datetime.datetime.now()


print('All Data Extracted\nTime Taken: %s\n\nCreating Figure'%(datetime.datetime.now()-startTime))
# Once all metrics are stored in a dictionary, we want to assemble arrays so that we can plot
# our data

# The first three arrays are for ages where multiple data points are associated with a particular
# age. The second two are for ages where there was only one data point associated with that age.
# The reason we differentiate between ages that have multiple data points for the metrics
# associated with them and those with only one is because we can't get the standard error if there's
# only one data point. Ages that have only one data point associated with them are plotted as light blue,
# slightly smaller dots than those where an average was possible.
x = []
y = []
e = []

x_1 = []
y_1 = []

# We sort through the dictionary that we compiled and add the data points to the arrays as we go
for age in AgeDictionary:
    if len(AgeDictionary[age]) != 1:
        standardError = stats.sem(AgeDictionary[age])
        if age != None and np.mean(AgeDictionary[age]) != None:
            x.append(age)
            y.append(np.mean(AgeDictionary[age]))
            e.append(standardError)
    else:
        x_1.append(age)
        y_1.append(AgeDictionary[age][0])

# Once all the arrays are complete, we plot them as a scatterplot
plt.errorbar(x, y, e, linestyle='None', markersize=5,marker='o')
plt.scatter(x_1, y_1, c = '#3FC0DA' ,linestyle = 'None', marker='.')
plt.xlabel('Age (MY)')
plt.ylabel('%s'%MetricsLabel)
plt.title('Full Proteins - %s'%MetricsLabel)
plt.savefig('Full Proteins_%s'%Metric.replace('-',''))
plt.close()
     
print("Done!\nTime Taken: %s"%(datetime.datetime.now()-startTime))

        
