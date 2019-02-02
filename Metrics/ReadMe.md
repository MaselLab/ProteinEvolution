
# Metrics
---------------------------------------------------------------------------------------

This repository is used to store all scripts designed to process and produce metrics for protein sequences and genomes. These metrics include (but are not limited to) predictions for: 
- Intrinsic structural disorder
- Aggregation propensity
- Clustering and dispersion
- Amino acid composition
- Codon Adaptation Index
 
 The scripts are listed below in alphabetical order.
 
-----------


## CalculateMeanPfamISDOverAllSpecies
-----

This script is used to calculate the mean ISD for each Pfam in the database PFAMphylostratigraphy. It uses two databases to do this, specifically:

1) EnsemblGenomes_Metrics_Complete
2) NCBIGenomes_Metrics_Complete

For each row in these datatables, the script pulls the Pfam UID and the domain ISD value associated with it and stored all values in a large dictionary. After every row has been processed, the script calculates mean and variance of each Pfam. These results are uploaded to the data table PfamUIDsTable_EnsemblAndNCBI. These values are not differentiated by species.


## RunIUPred2_IntrinsicDisorderPrediction
-----

  This script is used to generate intrinsic structural disorder predictions for proteins stored in a MySQL database. The results for each protein are uploaded to the same MySQL table they were extracted from
  
##### Dependencies
The user will need to run the script with Python3 with the following modules:

	 BioPython: https://anaconda.org/conda-forge/biopython
	 mysql.connector : https://pynative.com/install-mysql-connector-python/

To get Python3, anaconda can be downloaded from: https://www.anaconda.com/download/


The user will also need to have the IUPred2 executables. These can be accessed here by filling out a request form: https://iupred2a.elte.hu/download (Select IUPred2A to download)



