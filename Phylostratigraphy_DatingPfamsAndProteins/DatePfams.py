import os, sys, json, csv, copy, mysql.connector, datetime
from ete3 import Tree

'''
The purpose of this script is to assign dates to Pfams in our dataset.

First, in order to aquire the various dates associated with the divergence of all species in our dataset, we uploaded the list of species in our dataset to TimeTree.org. The Newick file was then downloaded.

Next, a database was made with each unique Pfam assigned to its own row. The databases containing all the genomes for this work were then searched for occurences of each Pfam. The list of all species in which the Pfam appeared was saved in the Pfam datatable.

Each Pfam 


General Age-Assigning Algorithm:
--------------------------------

The reason a different algorithm is used when only one species is present (algorithm [2]) is because the function [species_list].get_common_ancestor is used to find the subtree in the overall tree that has the common ancestor of all the species in the species list as the root node. If only one species is present, then the entire tree is returned as the subtree, invalidating the methodology used for algorithm [1] (half the length of the entire tree will always be assigned as the age of the Pfam in this case, so all Pfams with only one species associated with them will get the same, inaccurate age). 

-----

algorithm [1]: Used when more than one species is present:
   - The ete3 function [species_list].get_common_ancestor() is used to pull the subtree which has the common ancestor of all species in the species_list as the root node. 

   - The ete3 function subtree.iter_descendants("preorder") is used to iterate through the tree, starting at the root common node and going down one branch of the subtree. Because the Newick tree downloaded from TimeTree.org doesn't assign absolute ages to each node but rather relative ages, the length of each branch needs to be added to a total age variable. 

   - During the iteration process, it's checked to see whether the current position in the tree is a leaf or a node. If it's a node, we continue searching through the tree. As soon as a leaf is reached, the age of that leaf is added to the total age and we stop searching the tree. This is because we've completely gone down one branch of the subtree and have obtained the age. 

   - We then use the ete3 subtree.traverse("postorder") to find the age of the common ancestor in the tree. We divide this in half and add it to the total. This is the final age we upload to MySQL


algorithm [2]: Used when only one species is present

   - The algorithm for this case is much simpler than algorithm [1]. The species is located within the overall tree, the distance (age) is determined using the simple ete3 function node.dist and the age is divided in half. This is the age that is assigned to the Pfam
'''


##################################################
#      Input User-Specific Information Here      #
##################################################

start_time = datetime.datetime.now()
print('Script Executing\nCurrent Time: %s\n'%datetime.datetime.now())


# MySQL connection information
Database = 'PFAMphylostratigraphy'
User = ''
Host = ''
Password = ''

# Data table where Pfams are stored
PfamTable = 'PfamUIDsTable_EnsemblAndNCBI'

# Column names in Pfam table
PfamUIDColumn = 'PfamUID'
SpeciesColumn = 'Species'
FalseNegativeSpeciesColumn = 'Species_FalseNegatives'
AgeColumn = 'Age_MY'

# Newick tree used to date the pfams
NewickTreeFilename = 'PhylogeneticTree_AllSpecies.nwk'

# Paul has done a lot of work on searching all full genomes in all six reading frames
# to determine which species are false negatives in our dataset. These species are stored
# in the pfam table in the same format as the original species list but in a separate column.
# The user has the option to include these in the dating of each domain by setting the following
# option to True. To exclude them, set the option to False.
IncludeFalseNegativeSpecies = True


##################################################
#            Program Executes Below              #
##################################################

# First, the Newick tree is loaded
t = Tree(NewickTreeFilename,format=1)

# And a connection is made to the SQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)


PfamExtractionStatement = "SELECT "+','.join([PfamUIDColumn,SpeciesColumn,FalseNegativeSpeciesColumn])+" FROM " + PfamTable
mycursor.execute(PfamExtractionStatement)
pfamResults = mycursor.fetchall()

# Each pfam is then analyzed so an age can be assigned to it.
for pfam in pfamResults:
    
    UID = pfam[0]

    # The list of species each pfam shows up in is compiled depending on whether species flagged as false
    # negatives will be included in the dating process
    BaseSpeciesList = pfam[1].split(',')
    FalseNegativeSpeciesList = pfam[2].split(',') if pfam[2] != None else []
    speciesList = BaseSpeciesList + FalseNegativeSpeciesList if IncludeFalseNegativeSpecies == True else BaseSpeciesList
    numberOfSpecies = len(speciesList)
    
    # If more than one species exists for a particular pfam, then algorithm [1] is implemented (See: description in
    # the beginning of this script
    if numberOfSpecies != 1:
        # The subtree is pulled starting with the common ancestor of all species in the list as the root node
        subtree = t.get_common_ancestor(speciesList)
        # The total distance starts at zero and will be added to as each branch in the subtree is searched
        totalDistance = 0
        # We start at the base node and proceed down a branch of the subtree
        for node in subtree.iter_descendants("preorder"):
            # So long as our location is not a leaf, we add the relative age of the node to the total and continue to search the tree
            if node.is_leaf() == False:
                totalDistance += node.dist
            # Once we're located on a leaf, we have searched an entire branch of our subtree and, once we add the age of the leaf, we have acquired the
            # total age of our subtree (minus the age of the common ancestor)
            else:
                totalDistance += node.dist
                break
        # The variable maxLength is used to find the node that has the greatest number of children in the tree
        maxLength = 0
        # We then traverse the tree starting from the leaves and working our way toward the root node
        for node in subtree.traverse("postorder"):
            # We get the descendants of each node
            nodeDescendants = node.get_descendants()
            # We then compare the number of descendants to the greatest number found for a node thus far
            if len(nodeDescendants) > maxLength:
                # If it has more descendants than any thus far, we make a note of which node it was 
                rootNodeName = node.name
                # and find the age of that node
                rootNodeDistance = node.dist
        # The node with the greatest number of descendants in the common ancestor and we take its relative age and divide it in half (estimating the Pfam arose
        # halfway down that branch) and add it to the total
        totalDistance += (rootNodeDistance/2)
    else:
        # If only one species is associated with that Pfam, we implement algorithm [2]
        species = t&speciesList[0]
        totalDistance = (species.dist)/2
    # We then update the age of a Pfam UID in our MySQL table.
    AgeUpdateStatement = "UPDATE "+PfamTable+ " SET "+AgeColumn+"=%.4f WHERE "%totalDistance+PfamUIDColumn+" = '%s'"%UID
    mycursor.execute(AgeUpdateStatement)
    cnx.commit()
cnx.close()

print('Script Complete\nTime Taken: %s'%(datetime.datetime.now()-start_time))
sys.stdout.flush()

