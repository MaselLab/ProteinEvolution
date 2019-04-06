# Sequence Alignments

This repository stores any scripts that are used to create multiple sequence alignments

## Python Scripts

### CreatePfamAlignments_UploadToMySQL

Author : Sara Willis
Date   : March 8, 2019

This script is designed to read align pfams all sharing a common ID from many different species. It goes about this by accessing protein data tables (in our dataset there's one from ensembl and one from NCBI) and extracting protein sequences, the pfams that are contained within them protein sequence, the indices where they occur, and a quality check. Proteins that pass the quality check move forward and have the peptide sequences associated with their pfams extracted. 
There is the poosibility that a pfam will show up more than once in a given protein. If this is the case, the scripts looks at the lengths of the pfams and selects one at random from the upper 50th percentile of given lengths.
Once all the pfams meeting the previously described conditions, the script goes about the process of aligning all peptides associated with a particular pfam UID using Clustal Omega and uploads each aligned sequence to a MySQL data table. 
In some cases, Clustal Omega is unable to align all sequences associated with a particular pfam because there are too many. Thus far, it's not clear what to do with those. Another alignment program may be necessary. Until that's implemented, the script keeps track of these problematic pfams and prints them to a flat file for the user to look at once the script has been run.

DEPENDENCIES

This script requires Clustal Omega to run which can be [downloaded as a precompiled executable](http://www.clustal.org/omega/#Download).
The script also requires the non-standard python module mysql.connector

## Text Files

### PfamsTooLargeToAlign.csv

This csv file contains the IDs of all the Pfams which Clustal Omega was unable to align in the script CreatePfamAlignments_UploadToMySQL. They will need to be aligned using different methodology
