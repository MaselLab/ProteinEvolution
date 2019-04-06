# Phylostratigraphy_DatingPfamsAndProteins

## Python Scripts

### AssignProteinHomologyGroups_BasedOnOldestPfams

Author : Sara Willis
Date   : Thursday March 7, 2019

The purpose of this script is to assign full proteins to homology groups based on the oldest pfams that they are associated with. 
The way homology groups are assigned to proteins is to cluster proteins based on the PfamUID of the oldest pfam that shows up in that protein. The script would be relatively straight-forward if proteins and their pfams were well behaved and only one pfam with the oldest age showed up. Alas, it is not so. There are cases when there are two more more pfams that show up in a protein sequence that all share the same oldest pfam UID. As a consequence, we need a way to cluster pfams so homology groups can still be assigned to our proteins.
The clustering algorithm is as follows:

For simplicity, consider two pfam UIDs A and B. Suppose that they are the oldest pfams in their given protein. We then consider how frequently they show up with one another. If either:

      - Number of Occurences of A+B / Number of Occurences of A >= 50%
      - Number of Occurences of A+B / Number of Occurences of B >= 50%
    
A link is established. We then use single-link clustering to create pfam clusters. 
To generate the clusters, every pairwise comparison is made for all pfams that occur together and is output as a dictionary. This is saved as a flat file so this process doesn't necessarily need to be run each time the script is used. If a previous run has already generated the flat files, the script will just load previous results if the user sets the variable GeneratePfamClusters to False.
Once the pfam clusters are generated, the program uses the clustering information to assign homology groups to the full proteins in the database. It does this by using the following algorithm:

For simplicity, let's consider two clustering IDs. Lets say a protein has any number of oldest Pfams associated with it, we then look at the clustering IDs associated with those pfam IDs and assign homology groups based on those
If A and B are clustering IDs:

     - if A = B, then this is trivial and we assign the homology group number associated with that clustering ID
     - if A =/= B, then we generate a new homology group number associated with AB
   
We generalize for cases with n clustering ids. All homology IDs are stored in the protein metrics data tables.


### DateFullGenesUsingOldestPfam

Author       : Sara Willis
Date         : February 1, 2019
Updated      : April 5, 2019
Update Notes : Changed so MySQL column names defined only once allowing for greater ease-of-use

This script is used to date full proteins by assigning the age of the oldest Pfam that it contains. 
In order to run this script, the following databases are required:

         - Two databases containing protein sequences and the pfams that occur within them (one for NCBI, one for Ensembl)
         - Two databases where the ages of the protein sequences will be stores (one for NCBI, one for Ensembl)
         - A database where pfams are stored with their ages. The pfam UID should be the primary key in this database
         
The following databases are what have been used:

In the database PFAMphylostratigraphy, there are four data tables that contain full proteins and the metrics associated with them:

         1) EnsemblGenomes_Proteins_Complete
         2) EnsemblGenomes_ProteinMetrics_Complete
         3) NCBIGenomes_Proteins_Complete
         4) NCBIGenomes_ProteinMetrics_Complete
         
The Proteins data table contains a variety of data that include the UIDs of the pfams that show up in the body of the protein. The ProteinMetrics data tables contain the averages of the metrics that have been calculated for each protein as well as the age of the oldest Pfam that occurs in the protein. The UID of the ProteinMetrics data table is shared with the Proteins data table. 
For each of the source databases (Ensemble and NCBI), the PfamUIDs are extracted from the Proteins database for each protein, the oldest is found, and the age of the oldest is used to update the ProteinMetrics data table. 

In PFAMphylostratigraphy, the data table PfamUIDsTable_EnsemblAndNCBI contains all PfamUIDs and the ages that have been assigned to them.

### DatePfams

Author : Sara Willis
Date   : February 11, 2019

The purpose of this script is to assign dates to Pfams in our dataset.
First, in order to aquire the various dates associated with the divergence of all species in our dataset, we uploaded the list of species in our dataset to TimeTree.org. The Newick file was then downloaded.
Next, a database was made with each unique Pfam assigned to its own row. The databases containing all the genomes for this work were then searched for occurences of each Pfam. The list of all species in which the Pfam appeared was saved in the Pfam datatable.

General Age-Assigning Algorithm:

The reason a different algorithm is used when only one species is present (algorithm [2]) is because the function [species_list].get_common_ancestor is used to find the subtree in the overall tree that has the common ancestor of all the species in the species list as the root node. If only one species is present, then the entire tree is returned as the subtree, invalidating the methodology used for algorithm [1] (half the length of the entire tree will always be assigned as the age of the Pfam in this case, so all Pfams with only one species associated with them will get the same, inaccurate age). 

algorithm [1]: Used when more than one species is present:

- The ete3 function [species_list].get_common_ancestor() is used to pull the subtree which has the common ancestor of all species in the species_list as the root node. 
- The ete3 function subtree.iter_descendants("preorder") is used to iterate through the tree, starting at the root common node and going down one branch of the subtree. Because the Newick tree downloaded from TimeTree.org doesn't assign absolute ages to each node but rather relative ages, the length of each branch needs to be added to a total age variable. 
- During the iteration process, it's checked to see whether the current position in the tree is a leaf or a node. If it's a node, we continue searching through the tree. As soon as a leaf is reached, the age of that leaf is added to the total age and we stop searching the tree. This is because we've completely gone down one branch of the subtree and have obtained the age. 
- We then use the ete3 subtree.traverse("postorder") to find the age of the common ancestor in the tree. We divide this in half and add it to the total. This is the final age we upload to MySQL
   
algorithm [2]: Used when only one species is present

- The algorithm for this case is much simpler than algorithm [1]. The species is located within the overall tree, the distance (age) is determined using the simple ete3 function node.dist and the age is divided in half. This is the age that is assigned to the Pfam
         
## Newick Files

### PhylogeneticTree_AllSpecies.nwk
This is the phylogenetic tree of all species in our data set in Newick (machine readable) format 

## Text Files 

### SpeciesForPhylogeneticTree.txt
This is the list of species that has been used to generate the phylogenetic tree using [TimeTree](TimeTree.org)
