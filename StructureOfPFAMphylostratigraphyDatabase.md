# Structure of PFAMphylostratigraphy Database

PFAMphylostratigraphy is the MySQL database where ~400 fully-sequenced species' annotated genomes are stored. These genomes were collected for work on the Templeton grant from the Ensembl and NCBI databases subject to the following filters:

1. Species were listed as fully-sequenced in the [GOLD databases](https://gold.jgi.doe.gov/)
2. The genomes were fully-annotated so the coding sequences were available 
3. The species were listed in Timetree so ages could be assigned in a phylogenetic tree
4. Only genes were taken if they contained a Pfam domain 
  * In Ensembl, Pfam annotations were listed in Biomart. Ensembl uses InterProScan to assign Pfam annotations
  * The genes from NCBI had their Pfam annotations done locally using InterProScan so our assignment scheme is compatible 
5. Only genomes found in Ensembl that were present in BioMart were collected because these were the ones that had Pfam annotations

The scripts that were used in acquiring these genomes can be found in the GenomicDataCollectionScripts repository. 

The structure of this database is as follows:

## Species Table

There is a master list of all species collected. This is stored in the data table **SpeciesList** and contains the following data:
1. **SpeciesUID** -- This is the table-specific UID and is used throughout other tables to quickly link back to this table
2. **SpeciesName** -- The name of the species which can include the subspecies
3. **NewickSpeciesName** -- This is the name of the species in Newick format so that the species can be used in TimeTree and so dating can occur with any trees built in Newick format.
4. **SpeciesCommonName** -- This is only for some Ensembl species. This lists the common name of the species instead of the scientific name
5. **GoldStatus** -- This is the status of the species in the Gold Database, specifically if the sequence results are a permanent draft or are complete and published. All species collected were required to have at least "Draft" in their Gold entry and were not allowed to be listed as "Incomplete"
6. **StudyID** -- This is the ID of the study that sequenced the species
7. **GeneCount** -- The number of genes that appear in this specific species. *This needs to be updated*
8. **EnsemblAccession** -- Only for Ensembl species, NULL for NCBI
9. **SourceDatabase** -- This is the repository where the genomes were collected, can either be NCBI or one of five Ensembl databases including: 
  * Ensembl V93
  * Ensembl Metazoa
  * Ensembl Plants
  * Ensembl Fungi
  * Ensembl Protist
10. **GenebuildMethod** -- This is listed in Ensembl and is not available for NCBI
11. **InBiomart** -- This lists whether the species can be found in Biomart. If not, the genome is not included in our analyses
12. **OneToOne** -- Either 1 (True) or 0 (False). Signals whether there is a one-to-one correspondance between Protein IDs and Gene IDs, i.e. Genes don't have multiple transcripts (due to things like alternative splicing variants)
13. **Notes** -- In some cases, a genome was removed due to things like a species being a duplicate (accidentally collected from both Ensembl and NCBI, for example). These species were retained in the species list, but will have notes next to them explaining why they are no longer in the database. 
  


## Sequence Data, Protein Annotations, and Protein Metrics

There are two main groupings of tables that contain coding sequence, protein sequences, and metric data. One for NCBI, the other for Ensembl. They are named similarly for clarity. 

#### Coding Sequences

The coding sequences are stored in tables: 
* **EnsemblGenomes_Coding_Complete**
* **NCBIGenomes_Coding_Complete** 

These tables have almost identical layouts and contain the following entries:

1. **UID** -- table-specific UID
2. **Repository** -- This is the specific database the sequence was downloaded from. The databases are the following:
  * EnsemblV93
  * EnsemblMetazoa
  * EnsemblFungi
  * EnsemblMetazoa
  * NCBI
3. **SpeciesUID** -- This links the gene back to the SpeciesList table
4. **NewickSpecies** -- This is the species name as it appears in the master Newick-formatted tree file, excludes subspecies
5. **GeneID** -- This is the Gene ID listed in the source database
6. **TranscriptID** -- This is the specific transcript variant associated with the GeneID. There may be more than one associated with a particular gene. 
7. **ProteinID** -- Similar to the TranscriptID
8. [**TranscriptSupportLevel**](https://uswest.ensembl.org/Help/Glossary) -- This is only available for genes from Ensembl, and only for a very limited number so that it may be relatively useless. This is the experimental support for the existence of a particular transcript.
9. **CodingSequence** -- Protein coding genes only

### Protein Sequences
