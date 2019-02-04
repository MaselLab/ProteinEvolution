# Figures

This repository houses all scripts used to generate figures. 

## PlotAverageMetricVsAge_FullProteins.py

The purpose of this script is to plot the means and standard errors of various protein metrics vs. the age of that protein. The metrics determined for each protein are for that entire protein and not for any specific domains contained within it. Since there are a variety of metrics associated with proteins, this script has been designed so that the user can simply select which metric they wish to plot from the command line and the script will pull the relevant values. 

The metrics that can be selected from are the following 

   1) --ISD : Mean ISD predicted by IUPred 2
   2) --DensityOfAPRs : The number of aggregation-prone regions divided by the length of the protein calculated using the output of TANGO
   3) --DensityOfAAsInAPRs : The number of amino acids located in aggregation-prone regions divided by the length of the protein, calculated using the output of TANGO
   4) --Clustering_Trunc : The normalized index of dispersion using FILVM as the most hydrophobic amino acids and using one window
   5) --Clustering_AllFrames: The normalized index of dispersion using FILVM as the most hydrophobic amino acids, calculated using all possible windows
   6) --Length : The length of the amino acid sequence
   
#### TO RUN: 
        <python3 command> PlotAverageMetricVsAge_FullProteins.py --MetricOption

Choose metrics option from the above list or type --options to get a list in your terminal window
