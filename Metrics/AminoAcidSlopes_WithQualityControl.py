import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib
from scipy import stats
import numpy as np
matplotlib.use('agg')
import matplotlib.pyplot as plt

'''
The purpose of this file is to create a file on the changes in amino acid composition over time in species stored in MySQL databases

'''


##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################
Database = ''
User     = ''
Host     = ''
Password = ''

EnsemblDataTable = 'EnsemblGenomes_DomainMetrics_Complete'
NCBIDataTable = 'NCBIGenomes_DomainMetrics_Complete'


RunRCommand = '~/R-3.5.2/bin/Rscript'


Verbose = True


##########################################################################################
#                               Program Executes Below                                   #
##########################################################################################

startTime = datetime.datetime.now()
# Connects to MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

PfamDictionary = {}
AgeDictionary = {}
PfamAgesDictionary = {}

# Assembles a dictionary so each pfam UID points to its age
PullPfamAgesStatement = "SELECT PfamUID,Age_MY FROM PfamUIDsTable_EnsemblAndNCBI"
mycursor.execute(PullPfamAgesStatement)
PfamAges = mycursor.fetchall()
for PfamAge in PfamAges:
    PfamUID = PfamAge[0]
    Age = PfamAge[1]
    PfamAgesDictionary[PfamUID] = Age
    
# The order the amino acid fractions are stored in MySQL is
#      A,R,N,D,C,E,Q,G,H,O,I,L,K,M,F,P,U,S,T,W,Y,V
# So we define a list so that we can pull the correct values for each amino acid
AminoAcids = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

# The program will pull values from all relevant databases
DataTables = [EnsemblDataTable,NCBIDataTable]

for DataTable in DataTables:
    currentTime = datetime.datetime.now()
    if Verbose == True:
        print('Extracting from %s'%DataTable)
        sys.stdout.flush()

    currentTime = datetime.datetime.now()

    # Extracts the relevant metrics from the data table
    ExtractionStatement = "SELECT PfamUID,PercentAminoAcidComposition,PassedIUPredFilters FROM "+DataTable
    mycursor.execute(ExtractionStatement)
    results = mycursor.fetchall()
    # We only keep proteins that are of high quality (passed filters = True)
    for result in results:
        if result != None:
            PfamUID = result[0]
            PassedFilters = result[2]
            if PassedFilters != None:
                if result[1] != None and int(PassedFilters) != int(False):
                    AAComposition = [float(j) for j in result[1].split(',')]
                    # We then assemble a dictionary of the form:
                    # {Pfam : {AA_1 : [AA_1_comp_1,...,AA_1_comp_N],...,AA_M:[AA_M_comp_1,...,AA_M_comp_N]}}
                    # So, for example, if we have the Pfam Pfam0001, we'd find the fraction of each AA for each
                    # instance of that pfam. We then store each fraction in a list corresponding to each AA
                    # The overall fraction for that Pfam will then be the mean 
                    if PfamUID not in PfamDictionary:
                        PfamDictionary[PfamUID] = {}
                        for i in range(0,len(AminoAcids)):
                            PfamDictionary[PfamUID][AminoAcids[i]] = [AAComposition[i]]
                    else:
                        for i in range(0,len(AminoAcids)):
                            PfamDictionary[PfamUID][AminoAcids[i]].append(AAComposition[i])
    print('Extraction Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-currentTime))
    sys.stdout.flush()
    
print('All Extractions Complete\nTime Taken: %s\n\nAssembling Age Dictionary\n\n'%(datetime.datetime.now()-startTime))
sys.stdout.flush()
cnx.close()

currentTime = datetime.datetime.now()
# We then assemble create a dictionary where each pfam uid points to its age
for pfam in PfamDictionary:
    age = float(PfamAgesDictionary[pfam])
    if age not in AgeDictionary:
        AgeDictionary[age] = {}
        for i in range(0,len(AminoAcids)):
            AA = AminoAcids[i]
            AgeDictionary[age][AA] = [np.mean(PfamDictionary[pfam][AA])]
    else:
        for i in range(0,len(AminoAcids)):
            AA = AminoAcids[i]
            AgeDictionary[age][AA].append(np.mean(PfamDictionary[pfam][AA]))

print('Age Dictionary Created\nTime Taken: %s\n\nGenerating Linear Models\n\n'%(datetime.datetime.now()-currentTime))
sys.stdout.flush()

# We will save our results in a csv file
AminoAcidSlopesFile = open('AminoAcidSlopes.csv','w')

# Luke's table has the following format, so we save ours in the same format for the purposes of uniformity
AminoAcids = ['K','R','D','E','Q','S','N','T','P','Y','H','C','G','A','W','V','L','M','F','I']

# Each amino acid's results are turned into a data frame which is then box-cox transformed by R. A linear model
# is then fit to the model and the slope is saved
for AminoAcid in AminoAcids:
    DataFrameFile = open('DataFrame.csv','w')
    DataFrameFile.write('Age,AminoAcidFrequency\n')
    for age in AgeDictionary:
        if np.mean(AgeDictionary[age]['%s'%AminoAcid]) >0:
            DataFrameFile.write('%s,%s\n'%(age,np.mean(AgeDictionary[age]['%s'%AminoAcid])))
    DataFrameFile.close()

    RScriptFile = open('RunLinearModel.r','w')
    RScriptFile.write('library(MASS)\nlibrary(nlme)\nlibrary(lme4)\n')
    RScriptFile.write('df <- read.csv(file = "DataFrame.csv", header = T)\n')
    RScriptFile.write('bc <- boxcox(df$AminoAcidFrequency~1, lambda = seq(.1,.7,0.01))\n')
    RScriptFile.write('lambda <- bc$x[which.max(bc$y)]\n')
    RScriptFile.write('bc.transform <- function(x,L){x^L}\n')
    RScriptFile.write('AminoAcidFrequency.transform <- bc.transform(df$AminoAcidFrequency,lambda)\n')
    RScriptFile.write('df$AminoAcidFrequency.transform<-AminoAcidFrequency.transform\n')
    RScriptFile.write('SimpleLinearModel <- lm(df$AminoAcidFrequency.transform ~ df$Age)\n')
    RScriptFile.write('summary(SimpleLinearModel)\n')
    RScriptFile.close()
    os.system('%s ./RunLinearModel.r >ROutput.out'%RunRCommand)


    with open('ROutput.out','r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row != []:
                if row[0].split(' ')[0] == 'df$Age':
                    # The amino acid, its slope and the error are then written to the output data frame
                    AminoAcidSlopesFile.write('%s,%s,%s\n'%(AminoAcid,[i for i in row[0].split(' ') if i != ''][1],[i for i in row[0].split(' ') if i != ''][2]))

AminoAcidSlopesFile.close()
print('Slopes File Generated\n\nProgram Complete!Time Taken: %s'%(datetime.datetime.now()-startTime))
