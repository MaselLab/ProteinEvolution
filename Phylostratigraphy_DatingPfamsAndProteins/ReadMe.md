# Phylostratigraphy_DatingPfamsAndProteins

## AssignProteinHomologyGroups_BasedOnOldestPfams

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
