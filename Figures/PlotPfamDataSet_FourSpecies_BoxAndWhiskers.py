import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib
from scipy import stats
import numpy as np
matplotlib.use('agg')
import matplotlib.pyplot as plt


'''
Author : Sara Willis
Date   : February 1, 2019
-------------------------


TO RUN: 
        python3 <filename.py> --Metric

Choose metrics option from below or type --options


The purpose of this script is to plot the means and standard errors of various protein domain metrics vs. the age of that domain.

There are a variety of metrics that can be plotted. The metric is specified by the user as a command-line argument. They are the following:

   1) --ISD : The mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the amino acid sequence

'''


##########################################################################
#                               User Data                                #
##########################################################################

# MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

pathToRExecutable = "~/R-3.5.2/bin/Rscript"

ErrorBars = True
Loess = True

Verbose = True

# The SpeciesList UIDs for the species that will be plotted. Edit this to get a different sampling of species
SpeciesUIDs = [54,26,79,138]

# The MySQL database where Pfam UIDs are stored with their ages
PfamDatabase = "PfamUIDsTable_EnsemblAndNCBI"

NCBIDatabase = 'NCBIGenomes_DomainMetrics_Complete'
EnsemblDatabase = 'Genomes_Multicellular_DomainMetrics'


# If the user needs to make modifications to the names of the columns where the metrics are stored,
# or wishes to change the labeling scheme of their plot, they can make modifications to the dictionary
# below. The general format is:
#    Key = Input option
#    List entry 1: SQL column Name
#    List entry 2: Metric label on plot
MetricOptions = {'ISD': ['MeanISD_IUPred2_WithCys', 'Mean ISD','c(-0.000001,1)'],
                 'DensityOfAPRs': ['DensityOfAggProneRegions','Density of Aggregation-Prone Regions', 'c(-0.000001,.02)'],
                 'DensityOfAAsInAPRs': ['DensityOfAAsInAPRs',"Density of AAs in APRs",'c(-0.000001,.15)'],
                 'Clustering_Trunc': ['NormalizedIndexOfDispersion_Trunc_FILVM', 'Clustering','c(-0.000001,2)'],
                 'Clustering_AllFrames': ['NormalizedIndexOfDispersion_AllPhases_FILVM','Clustering (All Frames, FILVM)','c(-0.000001,2)'],
                 'Length': ['ProteinLength', 'Protein Length','c(-0.000001,2000)']}




##########################################################################
#                           Program Executes Below                       #
##########################################################################


# The program begins by trying to determine which metric the user would like plotted vs. time
# If the metric is located, the program proceeds. Otherwise, the program exits and notifies the
# user that there was an issue
try:
    Metric = sys.argv[1]
except:
    print('\n------------------\n\nInvalid format\n\nTo generated figures, a metric must be selected\n\nTo run, use the following format:\n\n' + '\33[34m' + 'python3 <filename.py> --Metric'+'\33[0m'+'\n\nFor help, use --options\n\n------------------\n')
    sys.exit(0)

if Metric == '--options':
    print('\n\nThe following options are available to plot vs. age:\n\n--ISD\n--DensityOfAPRs\n--DensityOfAAsInAPRs\n--Clustering_Trunc\n--Clustering_AllFrames\n--Length\n\n')
    sys.exit(0)
try:
    MetricsColumn = MetricOptions[Metric.replace('-','')][0]
    MetricsLabel = MetricOptions[Metric.replace('-','')][1]
    YLim = MetricOptions[Metric.replace('-','')][2]
except:
    print('\n\nInvalid metric option selected\n\nPlease specify a valid metric\n\nFor help, use --options\n\n')
    sys.exit(0)

startTime = datetime.datetime.now()

# The script then establishes a connection with the SQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We'll store each instance of a pfam in a dictionary, with each average metric value associated with
# that Pfam as a list. Once all the values associated with all the pfams are sorted, we'll take
# the average of the metrics for each Pfam
PfamDictionary = {}

# The program extracts all data values from the Ensembl and NCBI databases. 
SpeciesDatabases = [EnsemblDatabase,NCBIDatabase]
for SpeciesUID in SpeciesUIDs:
    print('\n\nExtracting Species %s\n\n'%SpeciesUID)
    for SpeciesDatabase in SpeciesDatabases:
        # We want to not only extract the metric, but also a quality control.
        # The quality control PassedIUPredFilters is an indication that a coding sequence
        # was low quality. If this is the case, we ignore the entry
        ExtractionStatement = "SELECT PfamUID,PassedIUPredFilters," + MetricsColumn + " FROM "+SpeciesDatabase + " WHERE SpeciesUID ="+str(SpeciesUID)
        print('\n\nExtracting data from %s\n\n'%SpeciesDatabase)
        mycursor.execute(ExtractionStatement)
        results = mycursor.fetchall()
        if results != None:
            for result in results:
                if result != None:
                    PfamUID,PassedFilters,ExtractedMetric = result
                    if ExtractedMetric != None and PassedFilters != int(False):
                        ExtractedMetric = float(ExtractedMetric)
                        if SpeciesUID not in PfamDictionary:
                            PfamDictionary[SpeciesUID] = {PfamUID:[ExtractedMetric]}
                        else:
                            if PfamUID not in PfamDictionary[SpeciesUID]:
                                PfamDictionary[SpeciesUID][PfamUID] = [ExtractedMetric]
                            else:
                                PfamDictionary[SpeciesUID][PfamUID].append(ExtractedMetric)
            print('Species Found and Extraction Complete\nTime Taken: %s\n\nCreating Age Dictionary\n\n'%(datetime.datetime.now()-startTime))
            currentTime = datetime.datetime.now()
            break

SpeciesNameDictionary = {}
for SpeciesUID in SpeciesUIDs:
    mycursor.execute("SELECT NewickSpeciesName FROM SpeciesList WHERE SpeciesUID ="+str(SpeciesUID))
    speciesName = mycursor.fetchone()
    NewickSpeciesName = speciesName[0]
    SpeciesNameDictionary[SpeciesUID] = NewickSpeciesName
# Each occurence of a metric is stored with all others of its age. The average value is found once all
# data entries are stored
AgeDictionary = {}

n = 0

for SpeciesUID in SpeciesUIDs:
    SpeciesName = SpeciesNameDictionary[SpeciesUID]
    TempDataFrameFilename = '%s.csv'%SpeciesName


    DataFrameForR = open(TempDataFrameFilename,'w')
    DataFrameForR.write('Metric,Age\n')

    for Pfam in PfamDictionary[SpeciesUID]:
        ExtractPfamAgeStatement = "SELECT Age_MY FROM "+PfamDatabase+ " WHERE PfamUID='%s'"%Pfam
        mycursor.execute(ExtractPfamAgeStatement)
        Age = mycursor.fetchone()[0]
        MeanPfamMetric = np.mean(PfamDictionary[SpeciesUID][Pfam])

        if Age not in AgeDictionary:
            AgeDictionary[Age] = [MeanPfamMetric]
        else:
            AgeDictionary[Age].append(MeanPfamMetric)
        if MeanPfamMetric == 0:
            MeanPfamMetric = float(1e-20)
        DataFrameForR.write('%s,%s\n'%(MeanPfamMetric,Age))
    DataFrameForR.close()

        

print('All Data Extracted\nTime Taken: %s\n\nCreating R Script'%(datetime.datetime.now()-startTime))

RScriptFilename = 'TempRScript.r'
RScript = open(RScriptFilename,'w')

RScript.write('library(MASS)\nlibrary(nlme)\nlibrary(lme4)\nlibrary(MuMIn)\n')
RScript.write('par(mfrow=c(2,2))\n')
for SpeciesUID in SpeciesUIDs:
    SpeciesName = SpeciesNameDictionary[SpeciesUID]
    TempDataFrameFilename = "%s.csv"%SpeciesName
    RScript.write('df_%s <- read.csv(file = "%s",header = T)\n'%(SpeciesName,TempDataFrameFilename))
    RScript.write('bc <- boxcox(df_%s$Metric~1, lambda = seq(.1,.7,0.01))\n'%SpeciesName)
    RScript.write('lambda <- bc$x[which.max(bc$y)]\n')
    RScript.write('bc.transform <- function(x,L){(x^L-1)/L}\n')
    RScript.write('bc.backtransform <- function(y,L){(L*y+1)^(1/L)}\n')
    RScript.write('Metric.transform <- bc.transform(df_%s$Metric,lambda)\n'%SpeciesName)
    RScript.write('df_%s$Metric.transform<-Metric.transform\n'%SpeciesName)


RScript.write('par(mfrow=c(2,2))\n')

for SpeciesUID in SpeciesUIDs:
    SpeciesName = SpeciesNameDictionary[SpeciesUID]
    
    RScript.write('SimpleLinearModel <- lm(df_%s$Metric ~ df_%s$Age)\n'%(SpeciesName,SpeciesName))
    RScript.write('TransformedLinearModel <- lm(df_%s$Metric.transform ~ df_%s$Age)\n'%(SpeciesName,SpeciesName))

    RScript.write('Intercept <- coef(SimpleLinearModel)["(Intercept)"]\n')
    RScript.write('Slope <- coef(SimpleLinearModel)["df_%s$Age"]\n'%SpeciesName)
    RScript.write('Pvalue <- summary(TransformedLinearModel)$coefficients[,4]  \n')
    RScript.write('BackTransformedIntercept <- bc.backtransform(Intercept,lambda)\n')
    RScript.write('attach(df_%s)\n'%SpeciesName)
    RScript.write('df_%s <- df_%s[order(Age),]\n'%(SpeciesName,SpeciesName))
    RScript.write('pvalueString <- paste(c("P=",signif(Pvalue,digits=2)[2]),sep = "",collapse="")\n')
    RScript.write('X <- df_%s$Age\nY <- df_%s$Metric\nTitle <- %s\n'%(SpeciesName,SpeciesName, '"'+SpeciesName.replace('_',' ') + '"'))
    RScript.write('loessMod10 <- loess(df_%s$Metric~df_%s$Age,span=1)\nsmoothed10 <- predict(loessMod10)\n'%(SpeciesName,SpeciesName))
    RScript.write('par(mar=c(5,4.5,3,1))\n')


    RScript.write('plot(df_%s$Age,df_%s$Metric,type="n",main=Title,xlab="Age (millions of years)",ylab="%s",cex.lab=1.5,col="dodgerblue",pch=20,cex.axis=1.4,cex.main=1.6,ylim=%s,xlim=c(-20,2300))\n'%(SpeciesName,SpeciesName,MetricsLabel,YLim))

    RScript.write('boxplot(Y~X,add=TRUE,at=c(unique(df_%s$Age)), boxfill = "lightblue",main="",xlab="",ylab="",xaxt="n",yaxt="n",position="dodge",pch=".",varwidth=TRUE,boxwex=20)\n'%SpeciesName)
    RScript.write('abline(Intercept,Slope,lwd=4)\n')
    RScript.write('usr <- par( "usr" )\n')
    RScript.write('pvalueString <- paste(c("p = ",signif(Pvalue,digits=2)[2]),sep = "",collapse="")\n')
    RScript.write('pvalueString\n')

    if Loess == True:
        RScript.write('lines(smoothed10,x=X,col="darkslategray4",lwd=4)\n')


    # Legend    
    #RScript.write('legend( usr[ 2 ]-(usr[2]*.7),usr[ 4 ]-(.05*usr[4]),legend=c("Linear Regression","Loess Regression"),col=c("black","darkslategray4"),lty=1:1,lwd=2,bg="white")\n')
    #RScript.write('myleg <- legend( usr[ 2 ]-(usr[2]*.65),usr[ 4 ]-(.05*usr[4]),legend=pvalueString,cex=1.4)\n')
    #RScript.write('rect(myleg$text$x[1]-diff(usr[1:2])/100 , myleg$rect$top-myleg$rect$h , 2295 , myleg$rect$top,col="white")\n') #myleg$rect$left+myleg$rect$w
    RScript.write('myleg <- legend( usr[ 2 ]-(usr[2]*.65),usr[ 4 ]-(.05*usr[4]),legend=pvalueString,cex=1.4)\n')
        
RScript.close()

TempROutputFilename = 'R.out'

os.system('%s %s'%(pathToRExecutable,RScriptFilename))

                



     
print("Done!\nTime Taken: %s"%(datetime.datetime.now()-startTime))

        
