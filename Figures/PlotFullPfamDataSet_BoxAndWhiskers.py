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

The output filename is rplot.pdf
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

# Setting Test to True will vastly reduce the number of sequences pulled to be plotted. This allows adjustments to be made to the plots
# without having to wait for the full database to download
Test = False

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
MetricOptions = {'ISD': ['MeanISD_IUPred2_WithCys', 'Mean ISD'],
                 'DensityOfAPRs': ['DensityOfAggProneRegions','Density of Aggregation-Prone Regions'],
                 'DensityOfAAsInAPRs': ['DensityOfAAsInAPRs',"Density of AAs in APRs"],
                 'Clustering_Trunc': ['NormalizedIndexOfDispersion_Trunc_FILVM', 'Clustering'],
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
    print('\n------------------\n\nInvalid format\n\nTo generated figures, a metric must be selected\n\nTo run, use the following format:\n\n' + '\33[34m' + 'python3 <filename.py> --Metric'+'\33[0m'+'\n\nFor help, use --options\n\n------------------\n')
    sys.exit(0)

if Metric == '--options':
    print('\n\nThe following options are available to plot vs. age:\n\n--ISD\n--DensityOfAPRs\n--DensityOfAAsInAPRs\n--Clustering_Trunc\n--Clustering_AllFrames\n--Length\n\n')
    sys.exit(0)
try:
    MetricsColumn = MetricOptions[Metric.replace('-','')][0]
    MetricsLabel = MetricOptions[Metric.replace('-','')][1]
    MetricChoice = copy.deepcopy(Metric.replace('-',''))
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

# The program extracts all data values from the Ensembl and NCBI databases
SpeciesDatabases = [EnsemblDatabase,NCBIDatabase]
for SpeciesDatabase in SpeciesDatabases:
    # We want to not only extract the metric, but also a quality control.
    # The quality control PassedIUPredFilters is an indication that a coding sequence
    # was low quality. If this is the case, we ignore the entry
    ExtractionStatement = "SELECT PfamUID,PassedIUPredFilters," + MetricsColumn + " FROM "+SpeciesDatabase #+ " WHERE UID < 100000"
    if Test == True:
        ExtractionStatement += " WHERE UID < 10000"
    print('\n\nExtracting data from %s\n\n'%SpeciesDatabase)
    mycursor.execute(ExtractionStatement)
    results = mycursor.fetchall()
    for result in results:
        if result != None:
            PfamUID,PassedFilters,ExtractedMetric = result
            if ExtractedMetric != None and PassedFilters != int(False):
                ExtractedMetric = float(ExtractedMetric)
                if PfamUID not in PfamDictionary:
                    PfamDictionary[PfamUID] = [ExtractedMetric]
                else:
                    PfamDictionary[PfamUID].append(ExtractedMetric)
    print('Extraction Complete\nTime Taken: %s\n\nCreating Age Dictionary\n\n'%(datetime.datetime.now()-startTime))
    currentTime = datetime.datetime.now()

# Each occurence of a metric is stored with all others of its age. The average value is found once all
# data entries are stored
AgeDictionary = {}

n = 0

TempDataFrameFilename = 'TempDataFrame.csv'


DataFrameForR = open(TempDataFrameFilename,'w')
DataFrameForR.write('Metric,Age\n')

for Pfam in PfamDictionary:
    ExtractPfamAgeStatement = "SELECT Age_MY FROM "+PfamDatabase+ " WHERE PfamUID='%s'"%Pfam
    mycursor.execute(ExtractPfamAgeStatement)
    Age = mycursor.fetchone()[0]
    MeanPfamMetric = np.mean(PfamDictionary[Pfam])
    
    if Age not in AgeDictionary:
        AgeDictionary[Age] = [MeanPfamMetric]
    else:
        AgeDictionary[Age].append(MeanPfamMetric)
    if MeanPfamMetric == 0:
        MeanPfamMetric = float(1e-20)
    DataFrameForR.write('%s,%s\n'%(MeanPfamMetric,Age))
DataFrameForR.close()

DataFrameForPlottingFilename = 'RDataFrameForPlotting'
DataFrameForPlotting = open(DataFrameForPlottingFilename,'w')
DataFrameForPlotting.write('Metric,Age\n')

for Age in AgeDictionary:
    for pfamValue in AgeDictionary[Age]:
        DataFrameForPlotting.write('%s,%s\n'%(pfamValue,Age))

DataFrameForPlotting.close()
        
        

print('All Data Extracted\nTime Taken: %s\n\nCreating R Script'%(datetime.datetime.now()-startTime))

RScriptFilename = 'TempRScript.r'
RScript = open(RScriptFilename,'w')

# Below, we're going to write the R script that will do the plotting for us

# We start by importing the packages we need
RScript.write('library(MASS)\nlibrary(nlme)\nlibrary(lme4)\nlibrary(MuMIn)\n')

# We read in the dataframe we've created
RScript.write('df <- read.csv(file = "%s",header = T)\n'%TempDataFrameFilename)
RScript.write('attach(df)\n')
RScript.write('df <- df[order(Age),]\n')
RScript.write('loessMod10 <- loess(df$Metric~df$Age,span=1)\nsmoothed10 <- predict(loessMod10)\n')

# The next job is to find the optimal lambda to transform our data with using a Box-Cox transform
RScript.write('bc <- boxcox(df$Metric~1, lambda = seq(.1,.7,0.01))\n')
RScript.write('lambda <- bc$x[which.max(bc$y)]\n')
RScript.write('bc.transform <- function(x,L){(x^L-1)/L}\n')
RScript.write('bc.backtransform <- function(y,L){(L*y+1)^(1/L)}\n')
RScript.write('Metric.transform <- bc.transform(df$Metric,lambda)\n')
RScript.write('df$Metric.transform<-Metric.transform\n')


    
# We then run a simple linear model on the untransformed metric to get the slope and intercept
RScript.write('SimpleLinearModel <- lm(df$Metric ~ df$Age)\n')
RScript.write('Intercept <- coef(SimpleLinearModel)["(Intercept)"]\n')
RScript.write('Slope <- coef(SimpleLinearModel)["df$Age"]\n')

# To get the pvalue, we run another linear model on the transformed data
RScript.write('LinearModel_transformed <- lm(df$Metric.transform ~ df$Age)\n')
RScript.write('Pvalue <- summary(LinearModel_transformed)$coefficients[,4]  \n')
RScript.write('pvalueString <- paste(c("P=",signif(Pvalue,digits=2)[2]),sep = "",collapse="")\n')

#RScript.write('df <-  read.csv(file = "%s",header = T)\n'%DataFrameForPlottingFilename)

# We let R know that we want to save the figure we create as a PDF with the filename in quotes
RScript.write('pdf("rplot.pdf",width=10)\n')
RScript.write('par(mar=c(5,4.5,3,1))\n') # Move the margins so we can see the axis labels
RScript.write('X <- df$Age\nY <- df$Metric\nTitle <- pvalueString\n')

# Generate a blank plot with axes labels, title, and margins so we can add everything on top of it (this is to get around some issues with mixing
# boxplots and linear regressions
if MetricChoice == 'ISD':
    RScript.write('plot(df$Age,df$Metric,type="n",main=Title,xlab="Age (millions of years)",ylab="%s",cex.lab=1.75,col="dodgerblue",pch=20,cex.axis=1.75,cex.main=1.7)\n'%MetricsLabel)
elif MetricChoice == 'Clustering_Trunc':
    RScript.write('plot(df$Age,df$Metric,type="n",main=Title,xlab="Age (millions of years)",ylab="%s",cex.lab=1.75,col="dodgerblue",pch=20,cex.axis=1.75,cex.main=1.7,ylim=c(-.0000001,2))\n'%MetricsLabel)
elif MetricChoice == 'DensityOfAPRs':
    RScript.write('plot(df$Age,df$Metric,type="n",main=Title,xlab="Age (millions of years)",ylab="%s",cex.lab=1.75,col="dodgerblue",pch=20,cex.axis=1.75,cex.main=1.7,ylim=c(-.0000001,0.02))\n'%MetricsLabel)

# Add the boxplot
RScript.write('boxplot(Y~X,add=TRUE,at=c(unique(df$Age)), boxfill = "lightblue",main="",xlab="",ylab="",xaxt="n",yaxt="n",position="dodge",pch=".",varwidth=TRUE,boxwex=20)\n')

# Add the linear regression
RScript.write('abline(Intercept,Slope,lwd=4)\n')



# If the user wants to plot a Loess Regression, we organize the data frame so it's possible,
# then determine the loess regression and plot it
if Loess == True:

    RScript.write('lines(smoothed10,x=X,col="darkslategray4",lwd=4)\n')
    
# Pulls the boundaries of the figure so we can standardize where the legend shows up
RScript.write('usr <- par( "usr" )\n')
RScript.write('legend( usr[ 2 ]-(usr[2]*.45),usr[ 4 ]-(.05*usr[4]),legend=c("Linear Regression","Loess Regression"),col=c("black","darkslategray4"),lty=1:1,cex=1.4,lwd=3)\n')

# Closes the PDF figure
RScript.write('dev.off()\n')
RScript.write('head(df)')
RScript.close()

TempROutputFilename = 'R.out'

# Runs the R script
os.system('%s %s'%(pathToRExecutable,RScriptFilename))
os.remove(RScriptFilename)

     
print("Done!\nTime Taken: %s"%(datetime.datetime.now()-startTime))

        
