# ReadMe: DataTables_CreationAndUpdates

## AssignIndicesToPfam.py

Author : Sara Willis
Date   : February 11, 2019

In each protein data table in the database PFAMphylostratigraphy, there is one row for for each protein and its associated Pfams. These tables have then been deconstructed to make the DomainMetrics data tables where there is one entry for each protein/pfam combo. What this script does is it gathers the UIDs in the DomainMetrics data tables, associates all the UIDs that belong to a particular protein, and uploads those back into the protein data table. This allows for unique keys that allow the user to go back and forth between data tables without ambiguity.



## CalculateMeanPfamISDOverAllSpecies.py

Created Friday January 4, 2019
Author: Sara Willis

The purpose of this script is to calculate the means and variances of all Pfam ISD values for all occurrences over all species in the PFAMphylostratigraphy database. The values are calculated by extracting from the metrics datatables:

    1) EnsembleGenomes_Metrics_Complete
    2) NCBIGenomes_Metrics_Complete
    
Once these values are calculated, they are uploaded into the data table PfamUIDsTable_EnsemblAndNCBI

## CreateOneEntryPerPfam_ProteinCombination

Author : Sara Willis
Date   : February 19, 2019

The purpose of this script is to take proteins and their associated Pfam domains from a protein data table and to split those entries so there is a separate entry for each protein/pfam combination. Those entries are then uploaded into a new data table

The data table that the entries are being loaded into should already exist prior to running the script.

Note: The UIDs in the Protein database should be integers

## CreatePfamTable_OneEntryPerPfam_ListAllAssociatedSpecies.py

This script is written to create a datatable that has a unique entry for each Pfam UID that appears in the Ensembl and NCBI datatables. 

It accomplishes this by extracting every domain entry in the datatables EnsemblGenomes_DomainMetrics_Complete and NCBIGenomes_DomainMetrics_Complete. In assembles these entries into a dictionary so that the metrics for each occurence of each PfamUID are stored in a list associated with that Pfam domain. 
