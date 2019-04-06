# Genomic Data Collection Scripts

## ExtractNCBIGenomes_UploadResultsToMySQL
Contains all scripts used to download coding sequences from NCBI. For more information, see the ReadMe in the repository's contents.

## DownloadMitochondrialAndChloroplastGeneIDs_Ensembl.py

Author : Sara Willis
Date   : February 11, 2019

This script is used to access the Ensembl databases so that it can extract IDs associated with chloroplast and mitochondrial genes.
This script was written so that mitochondrial and chloroplast genes could be identified in our dataset containing the coding sequences from over 400 genomes. This is in part because mitochondrial and chloroplast genes tend to use alternative coding alphabets. 
This script connects to the Ensembl FTP sites, checks various species names located in those databases and cross-references them with our databases. In the case that the species in the FTP site is one in our database, we extract that species, otherwise it's ignored. 
Once a species' Fasta file containing its coding sequences is downloaded, the description of each sequence is parsed. What we're looking for is either the identifier 'Pt' or 'MT' under 'chromosome:' that identifies the sequence as either chloroplast or mitochrondrial, respectively. 
If a chloroplast or mitochrondrial sequence is found, its UID is recorded in an output file for later use.


## DownloadMitochondrialAndChloroplastGeneIDs_NCBI.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this script is to search for mitochondrial or chloroplast genes in full genome fasta files on the NCBI FTP site. Only species in our dataset will be searched and all accession numbers associated with the chloroplast/mitochrondrial genes will be written to a file

