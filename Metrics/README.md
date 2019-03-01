
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


## CalculateCAI_AllSpecies
This script determines the [codon adaptation indices](https://en.wikipedia.org/wiki/Codon_Adaptation_Index) for all species stored in a MySQL database using their coding sequences and uploads them back into a user-specified data table in MySQL. 

Because every codon from every coding sequence (passing a quality filter) is used, and the fact that the CAI is calculated as a geometric mean of values associated with every codon, the products can become too small to be stored in Python, which forces values to 0 if they become too small. As a consequence, this script takes the log transform of the geometric mean to do all of its calculations before back-transforming, allowing us to use sums instead of products.

## CalculateDomainDisorder_IUPred.py

This script is intended to pull preexisiting raw IUpred scores from a full-protein run from a MySQL data table and to find the mean disorder prediction for each of its Pfam domains. The results get uploaded to a domain data table

To use this script, the IUPred scores should be stored in a database with comma-delimited strings of pfams IDs, their starting and stopping indices, and a list of UIDs associated with those pfam domains in a domain data table

**Note**: The UIDs for the tables in question should be integers, otherwise this script will exit with an error

## CalculateNormalizedIndexOfDispersion

The purpose of this script is to calculate the normalized indices of dispersion for protein sequences stored in a MySQL database. The functions used to calculate the index of dispersion were all created by Jason Bertram. The remainder of this script was written by Sara Willis.

## CalculatePercentAAComposition_FullProteins

This script is used to calculate the amino acid composition of protein sequences stored in a MySQL database. The script calculates the percent composition for each of the 20 amino acids, formats the values as a comma-delimited string, and updates a metrics table with the values. The metrics script should exist prior to running the script as this program will not create the data table for the user. 

**Note**: The UIDs associated with these tables should be integers, not strings. This is because, in general, it's a bit cost-intensive to pull gigantic databases all at once, if you're including all the protein sequences. This script skirts this problem by finding the max UID in the table and then iterating through so that only one entry is pulled at a time. If the UIDs are not integers, then this cannot be done.

The reason there are different scripts for pfam domains vs. full proteins is because there is not yet a table that contains all of the pfam sequences. If this is done, the domain-specific script can be deleted. Until then, it will be retained.

## CalculatePfamDomainAggregation_Tango.py


This script is intended to pull preexisiting raw Tango scores from a full-protein run from a MySQL data table and to find various aggregation metrics for the pfam domains associated with that protein sequence. The results will then be uploaded to a user-specified MySQL data table. 
To use this script, the Tango scores should be stored in a database with comma-delimited strings of pfams IDs, their starting and stopping indices, and a list of UIDs associated with those pfam domains in a domain data table

Pfam aggregation is calculated similarly to full protein aggregation. A Pfam is said to contain an aggregation-prone region if the boundaries of that pfam domain intersect with the boundaries of any distinct APR. 


**NOTE**: The UIDs for the tables in question should be integers, otherwise this script will exit with an error

## GCPercentContent_AllSpecies

The purpose of this script is to calculate the GC content for all species in a MySQL database using their coding sequences. The results are then uploaded to be stored with the relevant species

**Note**: This script initially calculated the GC content by dividing the number of G's + C's by the total number of nucleotides. It makes more sense to divide by the sum of the four nucleotides so we exclude unknown nucleic acids from our calculations. The script was updated on Github on March 1, 2019. 

**Note 2**: This script can be made more general by defining variables for the data tables in the user-options section. This may be done soon.



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



