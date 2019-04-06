
# Metrics
---------------------------------------------------------------------------------------

This repository is used to store all scripts designed to process and produce metrics for protein sequences and genomes. These metrics include (but are not limited to) predictions for: 
- Intrinsic structural disorder
- Aggregation propensity
- Clustering and dispersion
- Amino acid composition
- Codon Adaptation Index
 
 The scripts are listed in alphabetical order, sorted by programming language
 
 # Dependencies 
 
 Below are a list of the dependencies that may be required for the scripts included in this repository. 
### [Python3](https://www.anaconda.com/download/)
Necessary for a lot of scripts. This does not come preloaded on fusion, so user's of the lab will have to download their own local version

#### [BioPython](https://anaconda.org/conda-forge/biopython)
Perfect for translating sequences efficiently and reading sequence files (fasta, stockholm, etc...) as well as a host of other functions
 
#### [mysqlconnector](https://pynative.com/install-mysql-connector-python/)
A great python module that allows users to connect to a MySQL database 

### R

### Miscellaneous Executables

#### [IUPred2](https://iupred2a.elte.hu/download)

#### [Tango](http://tango.switchlab.org/)
Used for protein aggregation propensity

-----------

# Python Scripts

## CalculateCAI_AllSpecies

Author : Sara Willis
Date   : February 11, 2019

This script determines the [codon adaptation indices](https://en.wikipedia.org/wiki/Codon_Adaptation_Index) for all species stored in a MySQL database using their coding sequences and uploads them back into a user-specified data table in MySQL. 

Because every codon from every coding sequence (passing a quality filter) is used, and the fact that the CAI is calculated as a geometric mean of values associated with every codon, the products can become too small to be stored in Python, which forces values to 0 if they become too small. As a consequence, this script takes the log transform of the geometric mean to do all of its calculations before back-transforming, allowing us to use sums instead of products.

## CalculateDomainDisorder_IUPred.py

Author : Sara Willis
Date   : February 11, 2019

This script is intended to pull preexisiting raw IUpred scores from a full-protein run from a MySQL data table and to find the mean disorder prediction for each of its Pfam domains. The results get uploaded to a domain data table

To use this script, the IUPred scores should be stored in a database with comma-delimited strings of pfams IDs, their starting and stopping indices, and a list of UIDs associated with those pfam domains in a domain data table

**Note**: The UIDs for the tables in question should be integers, otherwise this script will exit with an error


## CalculatePercentAAComposition_FullProteins

Written : Friday January 11, 2019
Updated : Tuesday February 19, 2019
Author  : Sara Willis

This script is used to calculate the amino acid composition of protein sequences stored in a MySQL database. The script calculates the percent composition for each of the 20 amino acids, formats the values as a comma-delimited string, and updates a metrics table with the values. The metrics script should exist prior to running the script as this program will not create the data table for the user. 

**Note**: The UIDs associated with these tables should be integers, not strings. This is because, in general, it's a bit cost-intensive to pull gigantic databases all at once, if you're including all the protein sequences. This script skirts this problem by finding the max UID in the table and then iterating through so that only one entry is pulled at a time. If the UIDs are not integers, then this cannot be done.

The reason there are different scripts for pfam domains vs. full proteins is because there is not yet a table that contains all of the pfam sequences. If this is done, the domain-specific script can be deleted. Until then, it will be retained.

## CalculatePfamDomainAggregation_Tango

Author : Sara Willis
Date   : February 11, 2019

This script is intended to pull preexisiting raw Tango scores from a full-protein run from a MySQL data table and to find various aggregation metrics for the pfam domains associated with that protein sequence. The results will then be uploaded to a user-specified MySQL data table. 
To use this script, the Tango scores should be stored in a database with comma-delimited strings of pfams IDs, their starting and stopping indices, and a list of UIDs associated with those pfam domains in a domain data table

Pfam aggregation is calculated similarly to full protein aggregation. A Pfam is said to contain an aggregation-prone region if the boundaries of that pfam domain intersect with the boundaries of any distinct APR. 


**NOTE**: The UIDs for the tables in question should be integers, otherwise this script will exit with an error

## GCPercentContent_AllSpecies

Author   : Sara Willis
Date     : February 18, 2019
Modified : March 1, 2019 by Sara Willis

The purpose of this script is to calculate the GC content for all species in a MySQL database using their coding sequences. The results are then uploaded to be stored with the relevant species

**Note**: This script initially calculated the GC content by dividing the number of G's + C's by the total number of nucleotides. It makes more sense to divide by the sum of the four nucleotides so we exclude unknown nucleic acids from our calculations. The script was updated on Github on March 1, 2019. 

**Note 2**: This script can be made more general by defining variables for the data tables in the user-options section. This may be done soon.



## JasonsHydrophobicRunLengthFunctions

Author: Jason Bertram

These are a series of functions written by Jason Bertram to calculate the mean hydrophobic run lengths of amino acid sequences. Inside is also defined a hydrophobicity map -- a dictionary storing the scores associated with each of the amino acids telling the functions which amino acids are assumed to be the most hydrophobic. 


## RunClustering

The purpose of this script is to calculate the normalized indices of dispersion for protein sequences stored in a MySQL database. 
The functions used to calculate the index of dispersion were all created by Jason Bertram. The remainder of this script was written by Sara Willis.
EDITS
Non-Standard Amino Acid Abbreviations: There are some non-standard amino acid abbreviations that will cause the original function to fail. 
Specifically, B,O,U,J,Z. Sara has added these in to the original hydrophobicity dictionary so this script doesn't exit with an error.
Non-Standard Residue Abbreviations:
 
   B -- either D or N
   J -- either I or L
   Z -- either E or Q
   U -- Selenocysteine
   O -- Pyrrolysine
Note: Only the FILVM and FILVMW options in the dictionary have been updated. 
Removing Stop Codons: Stop codons, * , are not represented in the hydrophobicity dictionary and so need to be removed from any proteins being processed
to avoid exiting with a key error. A command has been added to this script that does this, but it should be kept in mind that this needs to be done
for any future scripts that make use of these functions
**IMPORTANT**
A note about speed: If you're updating a table and the protein sequence UID is not the primary key, make sure to convert the protein sequence UID to an index if the table is large. If the protein sequence UID is not indexed, the upload process could be *exceptionally* lengthy



## RunIUPred2

Author : Sara Willis
Date   : February 11, 2019

  This script is used to generate intrinsic structural disorder predictions for proteins stored in a MySQL database. The results for each protein are uploaded to the same MySQL table they were extracted from
  
##### Dependencies
The user will need to run the script with Python3 with the following modules:

	 BioPython: https://anaconda.org/conda-forge/biopython
	 mysql.connector : https://pynative.com/install-mysql-connector-python/

To get Python3, anaconda can be downloaded from: https://www.anaconda.com/download/


The user will also need to have the IUPred2 executables. These can be accessed here by filling out a request form: https://iupred2a.elte.hu/download (Select IUPred2A to download)


## RunIUPred2_StringOutputOnly

Author : Sara Willis
Date   : February 11, 2019

This is very similar to RunIUPred2, but instead only uploads the raw string to the destination table. These two scripts could be merged at some point to make the script RunIUPred2.py much more elegant and versatile.

## RunTangoInParallel

Author : Sara Willis
Date   : February 11, 2019

The purpose of this script is to read in protein sequences from a MySQL database and to calculate aggregation metrics for them using the executable Tango. 

This script makes use of the multiprocessing module to run Tango on multiple proteins in the database simultaneously. This allows the script to be run on fairly large datatables in a reasonable period of time. 

## ScrambleSequences_RunTangoInParallel

Author : Sara Willis
Date   : March 11, 2019

This script was written with the intention of being used to generated delta Tango scores. It reads in protein sequences from a MySQL database to generate randomly scrambled peptide sequences from those proteins and calculates aggregation metrics for the randomized sequences using the executable Tango. 

The scrambled peptides will be saved with their Tango output scores in a MySQL data table along with the UID that corresponds to the protein the randomized peptide was generated from. The user has the option to generate more than one randomized peptide sequence per protein using the NumberOfScrambledSequences variable in UserOptions. The table where the scrambled peptides and their metrics may or may not exist prior to running this script. If the user wants the script to generate the scrambled peptide table for them, they should set CreateRandomizedPeptideTable to True. A table with the name defined for the variable ScrambledTable in UserOptions will be created with the column names specified for that table (also in UserOptions).

This script makes use of the multiprocessing module to run Tango on multiple proteins in the database simultaneously. This allows the script to be run on fairly large datatables in a reasonable period of time. 

#### Dependencies

The user will need the executable Tango in order for this script to run: http://tango.crg.es/
The user will also need two non-standard Python modules that can be installed with the conda command:
   - BioPython       : https://anaconda.org/anaconda/biopython
   - mysql.connector : https://anaconda.org/anaconda/mysql-connector-python

# R Scripts

## calculate_normalized_clustering

This script is designed to calculate a normalized index of hydrophobic clustering, similar to Jason's Python script. It adjusts average run length by the expected run length from amino acid composition.
