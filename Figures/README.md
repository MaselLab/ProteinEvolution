# Figures

This repository houses all scripts used to generate figures. 

## PlotAverageMetricVsAge_FullProteins.py

#### TO RUN: 
```
python3 <filename.py> --Metric
```
Choose metrics option from the list below or type --options to get a list in your terminal window


The purpose of this script is to plot the means and standard errors of various protein metrics vs. the age of that protein. The metrics determined for each protein are for that entire protein and not for any specific domains contained within it. Since there are a variety of metrics associated with proteins, this script has been designed so that the user can simply select which metric they wish to plot from the command line and the script will pull the relevant values. 

The metrics that can be selected from are the following 

   1) --ISD : Mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the amino acid sequence
   



## SpeciesSpecific_PlotAverageMetricVsAge_PfamOnly.py

#### TO RUN:
```
python3 <filename.py> --Metric [--y_min=VALUE --y_max=VALUE]
```
The purpose of this script is to generate scatterplots of the means and standard errors of a user-specified metric vs. Age for pfam domains across individual species. 

The user should specify the metric that should be plotted in the command line using the format specified under TO RUN. The options are the following:

   1) --ISD : The mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the domain

The user also has the option to make the various plots have uniform axes by specifying a y_min and y_max value for each figure. These are optional arguments. If none are specified, matplotlib uses the default settings. 
