import os, sys, csv, mysql.connector, datetime
from scipy import stats
import numpy as np


'''
                                ==============================
				  Generate Amino Acid Slopes 
				==============================
Author : Sara Willis
Date   : March 14, 2019


This script is designed to generate amino acid slopes either for full genes or pfam domains.
This is done by collecting all proteins or domains, sorting them either by a homology group ID 
or by pfam ID, taking the average of all data points in a particular grouping, and creating a
data frame for all datapoints vs. their age. The data are either then normalized using a box-cox
transformation or not (this depends on the user input) and a linear regression is performed
on the data. This is done for each amino acid and the slopes and error terms are extracted from
each linear regression. All slopes and error terms are saved to a csv file for the user with the
filename prefix "AminoAcidSlopes" followed by user-specific options.

All data should be stored in a MySQL database for extraction.

This script allows the user to enter various options to customize the dataset they are interested in.
They may choose the kingdom, the transmembrane status, domains vs. proteins, etc. that they want to analyze.
All options are listed for the user below.

Note: order is not important for the command-line arguments
-----------------------------------------------------------



Usage            :  python3 ScriptName.py pfam|fullgene [test] [speciesuid=<n>] [kingdom=<x>] [transformed] [transmembrane=<x>]

pfam|fullgene    :  (required) pfam or fullgene. Specifies which dataset used to generate slopes                        
test             :  if included, only extracts first 1000 columns from each data table. For use with testing and debugging
speciesuid=n     :  Amino acid slopes for species with species UID n (found in SpeciesList MySQL table)                 
kingdom=x        :  Amino acid slopes for kingdom x (fungi, invertebrate, vertebrate, or plant)                         
transformed      :  Uses Box-Cox transform to normalize dataset before computing slope                                  
                    By default, if this argument is not included, the data will not be transformed                      
transmembrane=x  :  When not included as an argument, all data points are included.                                     
                    When x=true, only includes transmembrane proteins/domains                                           
                    When x=false, only includes proteins/domains predicted to not be transmembrane

                    
Domain Options: invertebrate, vertebrate, plant, fungi


A note about transforming the data:
===================================
Because Box-Cox transforms require values to be greater than zero, if the data are transformed, the value 0.5
is added to each data point. The reason this is required is because in some cases domains or proteins may be 
too short to include all amino acids and so sometimes the percent composition for a particular domain/protein
homology group is zero.



'''




##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################

# User's MySQL Connection Information
Database                      = ''
User                          = ''
Host                          = ''
Password                      = ''

# Full Protein Data Tables
EnsemblProteinDataTable       = 'Genomes_Multicellular_ProteinMetrics'
NCBIProteinDataTable          = 'NCBIGenomes_ProteinMetrics_Complete'

AgeColumn                     = 'AgeOfOldestPfam'
AACompositionColumn           = 'PercentAminoAcidComposition'
FilterColumn                  = 'PassedIUPredFilters'
HomologyGroupIDColumn         = 'HomologyGroupID'
TransmembraneColumn           = 'ExpAA'
UIDColumn                     = 'ProteinTableUID'
FullGeneSpeciesUIDColumn      = 'SpeciesUID'
ProteinLengthColumn           = 'ProteinLength'

# Pfam Data Tables
EnsemblPfamDataTable          = 'Genomes_Multicellular_DomainMetrics'
NCBIPfamDataTable             = 'NCBIGenomes_DomainMetrics_Complete'
PfamDataTableUIDColumn        = 'UID'
PfamDataTableAACompColumn     = 'PercentAminoAcidComposition'
PfamDataTableFilterColumn     = 'PassedIUPredFilters'
PfamDataTableProteinUIDColumn = 'ProteinTableUID'
PfamSpeciesUIDColumn          = 'SpeciesUID'
PfamDataTablePfamUIDColumn    = 'PfamUID'
DomainLengthColumn            = 'DomainLength'

# Pfam UID Table
PfamUIDTable                  = 'PfamUIDsTable_EnsemblAndNCBI'
PfamUIDColumn                 = 'PfamUID'
PfamAgeColumn                 = 'Age_MY'

# Species Table
SpeciesTable                  = 'SpeciesList'
SpeciesUIDColumn              = 'SpeciesUID'
KingdomColumn                 = 'Kingdom'

# User's path to their R executable 
RunRCommand                   = '~/R-3.5.2/bin/Rscript'

# How chatty the user wants the program to be 
Verbose = True


##########################################################################################
#                                     Submodules                                         #
##########################################################################################

'''
This submodule generates an R script to create a linear model for the generated data frame, extracts
the slopes and error terms and prints them to a csv file for the user
'''
def GenerateSlopes(HomologyDictionary,Transformed,filename):

    AminoAcidSlopesFile = open(filename,'w')
    AminoAcidSlopesFile.write('AA,Slope,Error\n')

    DataFrameFilename = 'DataFrame.csv'
    RScriptFilename = 'RunLinearModel.r'
    ROutputFilename = 'ROutput.out'


    # Luke's table has the following format, so we save ours in the same format for the purposes of uniformity
    AminoAcids = ['K','R','D','E','Q','S','N','T','P','Y','H','C','G','A','W','V','L','M','F','I']

    # Each amino acid's results are turned into a data frame which may or may not then be box-cox transformed by R, depending
    # on the user input arguments. A linear model  is then fit to the model and the slope and error are then saved to a csv file
    for AminoAcid in AminoAcids:
        print('Constructing Data Frame for %s'%AminoAcid)
        DataFrameFile = open(DataFrameFilename,'w')
        DataFrameFile.write('Age,AminoAcidFrequency\n')
        for entry in HomologyDictionary:
            # We find the mean AA frequency for each homology group (either full protein or pfam group, based on user options)
            MeanFrequency = np.mean(HomologyDictionary[entry][AminoAcid])
            Age = HomologyDictionary[entry]['Age']
            DataFrameFile.write('%s,%s\n'%(Age,MeanFrequency))
        DataFrameFile.close()

        # We then create an R script to run a linear model
        RScriptFile = open(RScriptFilename,'w')
        RScriptFile.write('library(MASS)\nlibrary(nlme)\nlibrary(lme4)\n')
        RScriptFile.write('df <- read.csv(file = "%s", header = T)\n'%DataFrameFilename)

        # We determine the optimal value of lambda using the shifted data values and perform a box-cox transform on the shifted
        # data
        if Transformed == True:
            RScriptFile.write('bc <- boxcox(df$AminoAcidFrequency~1, lambda = seq(.1,.7,0.01))\n')
            RScriptFile.write('lambda <- bc$x[which.max(bc$y)]\n')
            RScriptFile.write('bc.transform <- function(x,L){(x^L-1)/L}\n')
            RScriptFile.write('AminoAcidFrequency.transform <- bc.transform(df$AminoAcidFrequency,lambda)\n')
            RScriptFile.write('df$AminoAcidFrequency.transform<-AminoAcidFrequency.transform\n')
        # We make sure Age is treated as numeric so the linear model is fit appropriately
        RScriptFile.write('df$Age <- as.numeric(df$Age)\n')
        # The linear model is either fit to the transformed or untransformed data, determined by the user's input options
        if Transformed == True:
            RScriptFile.write('SimpleLinearModel <- lm(df$AminoAcidFrequency.transform ~ df$Age)\n')
        else:
            RScriptFile.write('SimpleLinearModel <- lm(df$AminoAcidFrequency ~ df$Age)\n')
        RScriptFile.write('summary(SimpleLinearModel)\n')
        RScriptFile.close()
        # We run the linear model and throw R's running output into a junk file so we don't clutter up our terminal. We save the rest of the output
        # to a file we then parse to get the slopes and error terms from
        os.system('%s ./%s &>%s'%(RunRCommand,RScriptFilename,ROutputFilename))
        with open(ROutputFilename,'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if row != []:
                    if row[0].split(' ')[0] == 'df$Age':
                        # The amino acid, its slope and the error are then written to the output data frame
                        print('Slope = %s\tError = %s\n'%([i for i in row[0].split(' ') if i != ''][1],[i for i in row[0].split(' ') if i != ''][2]))
                        AminoAcidSlopesFile.write('%s,%s,%s\n'%(AminoAcid,[i for i in row[0].split(' ') if i != ''][1],[i for i in row[0].split(' ') if i != ''][2]))
    # And we clean up the working directory
    FilesToDelete = [DataFrameFilename,RScriptFilename,ROutputFilename,'Rplots.pdf']
    for filename in FilesToDelete:
        if os.path.exists(filename) == True:
            os.remove(filename)
    AminoAcidSlopesFile.close()

# ---------------------------------------------------------------------------------------------------------------------------------

'''
This function is designed to create a dictionary of all full protein metrics, partitioned by homology group with subdictionaries designated
to the various amino acid frequencies associated with that homology group
'''

def CreateFullGeneDictionary(SpeciesUID = None,Test = False,Kingdom=None,Transmembrane=None,Transformed=False):
    
    HomologyDictionary = {}
    
    # The order the amino acid fractions are stored in MySQL is
    #      A,R,N,D,C,E,Q,G,H,O,I,L,K,M,F,P,U,S,T,W,Y,V
    # So we define a list so that we can pull the correct values for each amino acid
    AminoAcids = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

    # The program will pull values from all relevant databases. 
    DataTables = [EnsemblProteinDataTable,NCBIProteinDataTable]

    for DataTable in DataTables:
        currentTime = datetime.datetime.now()
        if Verbose == True:
            print('Extracting from %s'%DataTable)
            sys.stdout.flush()

        currentTime = datetime.datetime.now()


        # Extracts the relevant metrics from the data table
        ExtractionStatement = "SELECT "+','.join([AgeColumn,AACompositionColumn,FilterColumn,HomologyGroupIDColumn,TransmembraneColumn,ProteinLengthColumn])+" FROM "+DataTable

        # If a speciesUID is specified, we ignore the kingdom and only extract that species from whichever table contains it
        if SpeciesUID != None:
            ExtractionStatement += " WHERE "+FullGeneSpeciesUIDColumn+" = %s"%SpeciesUID
        # If a kingdom is specified and a speciesUID is not, then we find the UIDs of the species associated with that kingdom and extract the entries
        # that have those species UID values
        elif Kingdom != None:
            try:
                # We start by collecting the species UIDs assoicated with that kingdom from our species table
                KingdomExtractionStatement = "SELECT "+SpeciesUIDColumn+" FROM "+SpeciesTable+" WHERE "+KingdomColumn+"='%s'"%Kingdom
                mycursor.execute(KingdomExtractionStatement)
                KingdomUIDs = mycursor.fetchall()
                KingdomUIDs = [i[0] for i in KingdomUIDs]

            except:
                print('No kingdom with given specification exists. Check your syntax and try again\nFor help, use option -h')
                sys.exit(0)
            # We then add each species UID to our extraction statement with a where clause, separated by or statements so we get them all
            ExtractionStatement += " WHERE "
            ExtractionStatement += ' OR '.join([FullGeneSpeciesUIDColumn+'=%s'%i for i in KingdomUIDs])
            
        # If we're only testing the run, we only extract a subset of the entries from each table. This can't be used with species-specific or kingdom
        # runs since it's often the case that none of the values are extracted since the species UID may only have entries associated with higher UID values
        elif Test == True:
            ExtractionStatement += " WHERE "+UIDColumn+" < 1000"
        mycursor.execute(ExtractionStatement)
        results = mycursor.fetchall()
        # We only keep proteins that are of high quality (passed filters = True)
        for result in results:
            if result != None:
                Age = result[0]
                HomologyGroupID = result[3]
                PassedFilters = result[2]
                # ExpAA gives us a handle on whether a protein is transmembrane or not. Values > 18 are likely transmembrane, otherwise not
                ExpAA = float(result[4])
                Length = float(result[5])

                    
                if PassedFilters != None:
                    if result[1] != None and int(PassedFilters) != int(False):
                        # Transmembrane is always either True, False, or None and never changes allowing us to use the following statement to get
                        # uniform results
                        if (Transmembrane == True and ExpAA > 18) or (Transmembrane == False and ExpAA <= 18) or Transmembrane == None:
                            AAComposition = [float(j) for j in result[1].split(',')]
                            # We then assemble a dictionary of the form:
                            # {Pfam : {AA_1 : [AA_1_comp_1,...,AA_1_comp_N],...,AA_M:[AA_M_comp_1,...,AA_M_comp_N]}}
                            # So, for example, if we have the Pfam Pfam0001, we'd find the fraction of each AA for each
                            # instance of that pfam. We then store each fraction in a list corresponding to each AA
                            # The overall fraction for that Pfam will then be the mean 
                            if HomologyGroupID not in HomologyDictionary:
                                HomologyDictionary[HomologyGroupID] = {'Age':Age}
                                for i in range(0,len(AminoAcids)):
                                    if Transformed == True:
                                        AAComp = AAComposition[i]+(0.5/Length)
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]] = [AAComp]
                                    else:
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]] = [AAComposition[i]]
                            else:
                                for i in range(0,len(AminoAcids)):
                                    if Transformed == True:
                                        AAComp = AAComposition[i]+(0.5/Length)
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]].append(AAComp)
                                    else:  
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]].append(AAComposition[i])
        print('Extraction Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-currentTime))
        sys.stdout.flush()
    return HomologyDictionary


# ---------------------------------------------------------------------------------------------------------------------------------



def CreatePfamDomainDictionary(SpeciesUID = None,Test = False,Kingdom=None,Transmembrane=None,Transformed=False):
    if Transmembrane != None:
        DataTables = [EnsemblProteinDataTable,NCBIProteinDataTable]
        TransmembraneDictionary = {}
        
        for DataTable in DataTables:
            mycursor.execute("SELECT "+','.join([UIDColumn,TransmembraneColumn])+" FROM "+DataTable)
            TransmembraneResults = mycursor.fetchall()
            for TransmembraneResult in TransmembraneResults:
                ProteinTableUID, ExpAA = TransmembraneResult
                TransmembraneDictionary[ProteinTableUID] = float(ExpAA)
                
    AgeDictionary = {}
    mycursor.execute("SELECT "+','.join([PfamUIDColumn,PfamAgeColumn])+" FROM "+PfamUIDTable)
    PfamResults = mycursor.fetchall()
    for PfamResult in PfamResults:
        Pfam,Age = PfamResult
        AgeDictionary[Pfam] = Age
        
    
    HomologyDictionary = {}
    
    # The order the amino acid fractions are stored in MySQL is
    #      A,R,N,D,C,E,Q,G,H,O,I,L,K,M,F,P,U,S,T,W,Y,V
    # So we define a list so that we can pull the correct values for each amino acid
    AminoAcids = ['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']

    # The program will pull values from all relevant databases
    DataTables = [EnsemblPfamDataTable,NCBIPfamDataTable]

    for DataTable in DataTables:
        currentTime = datetime.datetime.now()
        if Verbose == True:
            print('Extracting from %s'%DataTable)
            sys.stdout.flush()

        currentTime = datetime.datetime.now()
    
        # Extracts the relevant metrics from the data table
        ExtractionStatement = "SELECT "+','.join([PfamDataTablePfamUIDColumn,PfamDataTableAACompColumn,PfamDataTableFilterColumn,PfamDataTableProteinUIDColumn,DomainLengthColumn])+" FROM "+DataTable
        if SpeciesUID != None:
            ExtractionStatement += " WHERE "+PfamSpeciesUIDColumn+" = %s"%SpeciesUID
        elif Kingdom != None:
            try:
                KingdomExtractionStatement = "SELECT "+SpeciesUIDColumn+" FROM "+SpeciesTable+" WHERE "+KingdomColumn+"='%s'"%Kingdom
                mycursor.execute(KingdomExtractionStatement)
                KingdomUIDs = mycursor.fetchall()
                KingdomUIDs = [i[0] for i in KingdomUIDs]

            except:
                print('No kingdom with given specification exists. Check your syntax and try again\nFor help, use option -h')
                sys.exit(0)
            ExtractionStatement += " WHERE "
            ExtractionStatement += ' OR '.join([PfamSpeciesUIDColumn+'=%s'%i for i in KingdomUIDs])

        elif Test == True:
            ExtractionStatement += " WHERE "+PfamDataTableUIDColumn+" < 1000"
        mycursor.execute(ExtractionStatement)
        results = mycursor.fetchall()
        # We only keep proteins that are of high quality (passed filters = True)
        for result in results:

            if result != None:
                Age = AgeDictionary[result[0]]
                HomologyGroupID = result[0]
                PassedFilters = result[2]
                ProteinTableUID = result[3]
                Length = float(result[4])
                if Transmembrane != None:
                    ExpAA = TransmembraneDictionary[ProteinTableUID]
                if PassedFilters != None:
                    if result[1] != None and int(PassedFilters) != int(False):
                        if (Transmembrane == True and ExpAA > 18) or (Transmembrane == False and ExpAA <= 18) or Transmembrane == None:
                            AAComposition = [float(j) for j in result[1].split(',')]
                            # We then assemble a dictionary of the form:
                            # {Pfam : {AA_1 : [AA_1_comp_1,...,AA_1_comp_N],...,AA_M:[AA_M_comp_1,...,AA_M_comp_N]}}
                            # So, for example, if we have the Pfam Pfam0001, we'd find the fraction of each AA for each
                            # instance of that pfam. We then store each fraction in a list corresponding to each AA
                            # The overall fraction for that Pfam will then be the mean 
                            if HomologyGroupID not in HomologyDictionary:
                                HomologyDictionary[HomologyGroupID] = {'Age':Age}
                                for i in range(0,len(AminoAcids)):
                                    if Transformed == True:
                                        AAComp = AAComposition[i]+(0.5/Length)
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]] = [AAComp]
                                    else:
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]] = [AAComposition[i]]
                            else:
                                for i in range(0,len(AminoAcids)):
                                    if Transformed == True:
                                        AAComp = AAComposition[i]+(0.5/Length)
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]].append(AAComp)
                                    else:
                                        HomologyDictionary[HomologyGroupID][AminoAcids[i]].append(AAComposition[i])
        print('Extraction Complete\nTime Taken: %s\n\n'%(datetime.datetime.now()-currentTime))
        sys.stdout.flush()
    return HomologyDictionary


# ---------------------------------------------------------------------------------------------------------------------------------

'''
This function will display to the user how to run the program and anything else about this scripts functionality. To access it, use the command -h on the command line
'''

def About():
    print('\n\n\t\t\t\t===============================\n\t\t\t\tGenerate Amino Acid Slopes File\n\t\t\t\t===============================\n\nThis script is designed to generate amino acid slopes either for full genes or pfam domains.\nThis is done by collecting all proteins or domains, sorting them either by a homology group ID \nor by pfam ID, taking the average of all data points in a particular grouping, and creating a\ndata frame for all datapoints vs. their age. The data are either then normalized using a box-cox\ntransformation or not (this depends on the user input) and a linear regression is performed\non the data. This is done for each amino acid and the slopes and error terms are extracted from\neach linear regression. All slopes and error terms are saved to a csv file for the user with the\nfilename prefix "AminoAcidSlopes" followed by user-specific options.\n\nThis script allows the user to enter various options to customize the dataset they are interested in.\nThey may choose the kingdom, the transmembrane status, domains vs. proteins, etc. that they want to analyze.\nAll options are listed for the user below.\n\nNote: order is not important for the command-line arguments')
    print('-----------------------------------------------------------\n\n\n')
    print('{:<17}{:<3}{:<100}'.format('Usage',':', 'python3 ScriptName.py pfam|fullgene [test] [speciesuid=<n>] [kingdom=<x>] [transformed] [transmembrane=<x>]')+'\n')
    print('{:<17}{:<3}{:<100}'.format('pfam|fullgene',':','(required) pfam or fullgene. Specifies which dataset used to generate slopes'))
    print('{:<17}{:<3}{:<100}'.format('test',':','if included, only extracts first 1000 columns from each data table. For use with testing and debugging'))
    print('{:<17}{:<3}{:<100}'.format('speciesuid=n',':','Amino acid slopes for species with species UID n (found in SpeciesList MySQL table)'))
    print('{:<17}{:<3}{:<100}'.format('kingdom=x',':','Amino acid slopes for kingdom x (fungi, invertebrate, vertebrate, or plant)'))
    print('{:<17}{:<3}{:<100}'.format('transformed',':','Uses Box-Cox transform to normalize dataset before computing slope'))
    print('{:<17}{:<3}{:<100}'.format('','','By default, if this argument is not included, the data will not be transformed'))
    print('{:<17}{:<3}{:<100}'.format('transmembrane=x',':','When not included as an argument, all data points are included.'))
    print('{:<17}{:<3}{:<100}'.format('','','When x=true, only includes transmembrane proteins/domains'))
    print('{:<17}{:<3}{:<100}'.format('','','When x=false, only includes proteins/domains predicted to not be transmembrane\n\n'))
    print('Domain Options: invertebrate, vertebrate, plant, fungi\n\n')
    print('A note about transforming the data:\n===================================\nBecause Box-Cox transforms require values to be greater than zero, if the data are transformed, the value 0.5\nis added to each data point. The reason this is required is because in some cases domains or proteins may be \ntoo short to include all amino acids and so sometimes the percent composition for a particular domain/protein\nhomology group is zero.\n\n\n')
    sys.exit(0)

    
# ---------------------------------------------------------------------------------------------------------------------------------

'''
A small function telling the user something has gone awry
'''

def ErrorInput():
    print('\n\nSyntax Error. A valid option must be selected.\nTo see user options, use the command-line argument -h\n\n')
    sys.exit(0)


# ---------------------------------------------------------------------------------------------------------------------------------


'''
This simple submodules gets all input user options in order for the rest of the program to proceed while configuring the output
filename
'''
def ParseUserOptions(arguments,filename_prefix):
    try:
        SpeciesUID = arguments['speciesuid']
        filename_prefix += '_SpeciesUID%s'%SpeciesUID
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[36;1mSpecies UID Selected : %s'%SpeciesUID,'--','This will override any kingdoms selected'))
    except:
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[33;1mNo species UID selected','--','Ignoring specific species'))
        SpeciesUID = None
    try:
        Transformed = arguments['transformed']
        filename_prefix += '_BoxCoxTransformed'
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[36;1mTransformation Selected','--','Linear models generated with Box-Cox transformed data'))
    except:
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[33;1mNo transformation selected','--','Linear models generated by untransformed data'))
        Transformed = False
        filename_prefix += '_NoTransform'
    try:
        Kingdom = arguments['kingdom']
        filename_prefix += '_%s'%Kingdom
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[36;1mKingdom selected','--',Kingdom))
    except:
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[33;1mNo kingdom selected','--','Including all kingdoms'))
        Kingdom = None

    try:
        Transmembrane = arguments['transmembrane']
        if Transmembrane == 'true':
            print('\n{:<45}{:<5}{:<100}'.format('\u001b[36;1mTransmembrane option selected','--','Only transmembrane proteins included'))
            Transmembrane = True
            filename_prefix += '_TransmembraneOnly'
        else:
            print('\n{:<45}{:<5}{:<100}'.format('\u001b[36;1mTransmembrane option selected','--','Only non-transmembrane proteins included'))
            filename_prefix += '_NonTransmembrane'
            Transmembrane = False
    except:
        print('\n{:<45}{:<5}{:<100}'.format('\u001b[33;1mNo transmembrane options specified','--','Including entire dataset'))
        Transmembrane = None
    try:
        Test = arguments['test']
        Test = True
    except:
        Test = False
    print('\u001b[0m\n\n')
    sys.exit(0)
    return SpeciesUID,Transformed,Kingdom,Transmembrane,Test,filename_prefix




##########################################################################################
#                               Program Executes Below                                   #
##########################################################################################


startTime = datetime.datetime.now()


# The program reads in all the user's arguments and assembles a dictionary that will direct the script how to run
arguments = {}
for i in sys.argv[1:]:
    if '=' in i:
        key,argument = i.split('=')
        arguments[key.replace('-','').lower()]=argument.lower()
    else:
        arguments[i.replace('-','').lower()] = True
        
# If the user has asked for help, the About function will run and the program will exit
if 'h' in arguments or 'help' in arguments or 'about' in arguments or 'options' in arguments:
    About()

# The script establishes a connection to MySQL
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# All amino acid output files start with the same prefix and have descriptors appended for each user option to
# differentiate files
filename_prefix = 'AminoAcidSlopes'

# We set up the variable names that will tell the program how to run
SpeciesUID,Transformed,Kingdom,Transmembrane,Test,filename_prefix = ParseUserOptions(arguments,filename_prefix)

# Pfams vs. full genes require slightly different functions and access to different databases, so we run them differently getting
# back dictionaries in the same format
if 'fullgene' in arguments:
    HomologyDictionary = CreateFullGeneDictionary(SpeciesUID,Test,Kingdom,Transmembrane,Transformed)
    filename_prefix += '_FullGenes'
elif 'pfam' in arguments:
    HomologyDictionary = CreatePfamDomainDictionary(SpeciesUID,Test,Kingdom,Transmembrane,Transformed)
    filename_prefix += '_PfamOnly'
# This script requires either pfam or fullgene to be specified to run, otherwise it quits with a complaint
else:
    ErrorInput()
    
filename = filename_prefix + '.csv'

# Once the dictionary is assembled, we generate linear models for each of the amino acids
GenerateSlopes(HomologyDictionary,Transformed,filename)


print('Amino Acid Slopes File Generated!\nTime Taken: %s'%(datetime.datetime.now()-startTime))


