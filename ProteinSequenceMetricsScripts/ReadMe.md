README

# ProteinSequenceMetricsScripts
---------------------------------------------------------------------------------------

This repository is used to store all scripts designed to process and produce metrics for protein sequences. These metrics include (but are not limited to) predictions for: 
- Intrinsic structural disorder
- Aggregation propensity
- Clustering and dispersion
- Amino acid composition
 
 The scripts are listed below in alphabetical order. All scripts are stored in their own subdirectory with a readme for ease-of-use.
 -----------

## RunIUPred2 
-----

  This script is used to generate intrinsic structural disorder predictions for proteins stored in a MySQL database. The results for each protein are uploaded to the same MySQL table they were extracted from
  
##### Dependencies
The user will need to run the script with Python3 with the following modules:

	 BioPython: https://anaconda.org/conda-forge/biopython
	 mysql.connector : https://pynative.com/install-mysql-connector-python/

To get Python3, anaconda can be downloaded from: https://www.anaconda.com/download/


The user will also need to have the IUPred2 executables. These can be accessed here by filling out a request form: https://iupred2a.elte.hu/download (Select IUPred2A to download)



