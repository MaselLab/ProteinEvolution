import os, json, sys, mysql.connector, datetime

'''

Author : Sara Willis
Date   : Thursday March 7, 2019

The purpose of this script is to assign full proteins to homology groups based on the oldest pfams that they are associated with. 

The way homology groups are assigned to proteins is to cluster proteins based on the PfamUID of the oldest pfam that shows up in that protein. The script would be relatively straight-forward if proteins and their pfams were well behaved and only one pfam with the oldest age showed up. Alas, it is not so. There are cases when there are two more more pfams that show up in a protein sequence that all share the same oldest age. As a consequence, we need a way to cluster pfams so homology groups can still be assigned to our proteins.

The clustering algorithm is as follows:

For simplicity, consider two pfam UIDs A and B. Suppose that they are the oldest pfams in their given protein. We then consider how frequently they show up with one another. If either:

    - Number of Occurences of A+B / Number of Occurences of A >= 50%
    or
    - Number of Occurences of A+B / Number of Occurences of B >= 50%

A link is established. We then use single-link clustering to create pfam clusters. 

To generate the clusters, every pairwise comparison is made for all pfams that occur together and is output as a dictionary. This is saved as a flat file so this process doesn't necessarily need to be run each time the script is used. If a previous run has already generated the flat files, the script will just load previous results if the user sets the variable GeneratePfamClusters to False.

Once the pfam clusters are generated, the program uses the clustering information to assign homology groups to the full proteins in the database. It does this by using the following algorithm:

For simplicity, let's consider two clustering IDs. Lets say a protein has any number of oldest Pfams associated with it, we then look at the clustering IDs associated with those pfam IDs and assign homology groups based on those
If A and B are clustering IDs:

   - if A = B, then this is trivial and we assign the homology group number associated with that clustering ID
   - if A =/= B, then we generate a new homology group number associated with AB 

We generalize for cases with n clustering ids.

All clustering IDs are stored in the protein metrics data tables.

'''

########################################################################
#                        User-Specific Information                     #                     
########################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# Pfam data table with ages. This is used to determine which are the oldest Pfams
# in a protein
PfamDataTable = 'PfamUIDsTable_EnsemblAndNCBI'
PfamUIDColumn = 'PfamUID'
AgeColumn = 'Age_Oldest_MY'

# The protein data tables where the pfam UIDs, metrics, and sequences are stored
EnsemblProteinTable = 'Genomes_Multicellular_Protein'
EnsemblProteinMetricsTable = 'Genomes_Multicellular_ProteinMetrics'
Ensembl_UIDColumn = 'UID'
Ensembl_PfamColumn = 'PfamUID'
Ensembl_HomGpColumn = 'HomologyGroupID_LUCA_LECA_Updated'
Ensembl_ProteinMetricsUIDColumn = 'ProteinTableUID'


NCBIProteinTable = 'NCBIGenomes_Protein_Complete'
NCBIProteinMetricsTable = 'NCBIGenomes_ProteinMetrics_Complete'
NCBI_UIDColumn = 'UID'
NCBI_PfamColumn = 'PfamUID'
NCBI_HomGpColumn = 'HomologyGroupID_LUCA_LECA_Updated'
NCBI_ProteinMetricsUIDColumn = 'ProteinTableUID'

# The user has the option to generate the Pfam clusters with this script.
# They must exist for the program to run in its entirety, so it needs to be executed
# at least once. Once it's run, however, the dictionaries that are generated by the run are dumped
# to flat files so this component doesn't need to be run every time. This is beneficial since the
# clustering process is somewhat time-intensive.
GeneratePfamClusters = True

# Below are the filenames containing the generated pfam clustering dictionaries
# These will be generated if GeneratePfamClusters is set to True, and will be loaded
# if it's set to False
Pfam_To_ClusteringID_Dictionary = 'PfamToClusteringID.dict'
ClusteringID_To_PfamCluster_Dictionary = 'ClusteringIDToPfam.dict'

Verbose = True



########################################################################
#                                Submodules                            #                     
########################################################################

'''
The purpose of this submodule is to cluster Pfams into groups based on the number of times they occur with one another

'''

def GeneratePfamClusteringDictionaries(PfamAgesDictionary,DataTablesDictionary,Verbose,PfamToClusteringID_Dictionary,ClusteringID_To_PfamCluster_Dictionary):
    startTime = datetime.datetime.now()
    # Establish a connection to the server
    cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
    mycursor = cnx.cursor(buffered = True)
    
    # We want to look at the proteins in each of our databases and we want to keep track of which is the oldest pfam
    # associated with that protein

    # Creates a list of sets, each set being a grouping of the oldest pfams that show up in a particular protein
    OldestPfams = []
    
    # PfamSet will be used to keep track of all the unique pfam UIDs that show up in our analyses
    PfamSet = set([])

    # SingleOccurences keeps track of the number of times a given pfam shows up alone
    SingleOccurences = {}

    # All data tables are searched
    for DataTableSource in DataTablesDictionary:
        if Verbose == True:
            print('Extracting %s Sequences'%DataTableSource)
            sys.stdout.flush()
        ProteinTable = DataTablesDictionary[DataTableSource]['ProteinTable']
        ProteinMetricsTable = DataTablesDictionary[DataTableSource]['ProteinMetricsTable']
        UIDColumn = DataTablesDictionary[DataTableSource]['UIDColumn']
        PfamColumn = DataTablesDictionary[DataTableSource]['PfamColumn']
        HomGpColumn = DataTablesDictionary[DataTableSource]['HomGpColumn']
        ProteinMetricsUIDColumn = DataTablesDictionary[DataTableSource]['ProteinMetricsUIDColumn']

        mycursor.execute("SELECT "+','.join([UIDColumn,PfamColumn]) + " FROM " + ProteinTable)
        if Verbose == True:
            print('Parsing %s Results\nTime Taken: %s\n'%(DataTableSource,datetime.datetime.now()-startTime))
            intermediateTime = datetime.datetime.now()
            sys.stdout.flush()
        ProteinResults = mycursor.fetchall()
        for proteinResult in ProteinResults:
            # This dictionary will be used to keep track of the oldest pfams for
            # a particular protein and so will be reset for each entry
            proteinSpecificPfams = {}
    
            UID, Pfams = proteinResult
            PfamList = Pfams.split(',')

            for Pfam in PfamList:
                if Pfam in PfamAgesDictionary:
                    AgeOfPfam = PfamAgesDictionary[Pfam]
                    if AgeOfPfam not in proteinSpecificPfams:
                        # We make sure we only have one copy for each pfam that shows up
                        # in the protein by making the collection of PfamUIDs a set
                        proteinSpecificPfams[AgeOfPfam] = set([Pfam])
                    else:
                        proteinSpecificPfams[AgeOfPfam].add(Pfam)
            # We then select the oldest age that shows up in a protein, assuming one exists
            if len(proteinSpecificPfams) != 0:
                MaxAge = max(proteinSpecificPfams)
                # If there are multiple pfams sharing that oldest age, we add the set of those Pfams
                # that share that particular age to the list OldestPfams. This will be used later to
                # determine the frequency each pfam shows up with others
                if len(proteinSpecificPfams[MaxAge]) != 1:
                    OldestPfams.append(proteinSpecificPfams[MaxAge])
                    PfamSet.update(proteinSpecificPfams[MaxAge])
                else:
                    MaxPfam = list(proteinSpecificPfams[MaxAge])[0]
                    PfamSet.add(MaxPfam)
                    if MaxPfam not in SingleOccurences:
                        SingleOccurences[MaxPfam] = 1
                    else:
                        SingleOccurences[MaxPfam] += 1

        if Verbose == True:
            print('%s Results Parsed\nTime Taken: %s'%(DataTableSource, datetime.datetime.now()-intermediateTime))
            sys.stdout.flush()
    
    # We will now construct a dictionary that keeps track of every pfam, the number of times
    # it shows up alone and the number of times it shows up with at least one other pfam sharing
    # a common age, as well as every pfam that shows up with it and the number of times that it
    # shows up with that particular pfam
    
    PfamGroupingDictionary = {}
    for pfam in PfamSet:
        PfamGroupingDictionary[pfam] = {}
        totalCombined = 0
        for pfamCluster in OldestPfams:
            if pfam in pfamCluster:
                totalCombined += 1
                for linkedPfam in pfamCluster:
                    if linkedPfam != pfam:
                        if linkedPfam not in PfamGroupingDictionary[pfam]:
                            PfamGroupingDictionary[pfam][linkedPfam] = 1
                        else:
                            PfamGroupingDictionary[pfam][linkedPfam] += 1
            PfamGroupingDictionary[pfam]['Total Combined']= totalCombined
            try:
                PfamGroupingDictionary[pfam]['Total Alone'] = SingleOccurences[pfam]
            except KeyError:
                PfamGroupingDictionary[pfam]['Total Alone'] = 0

    if Verbose == True:
        print('Pfam Grouping Dictionary Constructed\nTime Taken: %s\n\nGenerating Pfam Clusters'%(datetime.datetime.now()-intermediateTime))
        intermediateTime = datetime.datetime.now()
        sys.stdout.flush()
    # Now our goal is to sort pfams into groups which we achieve by single-link clustering.
    # To create a link between two pfams A and B, we require that at least one of the following conditions is met:
    #   1) Number of occurences of A+B / Number of occurences of A >= 50%
    #   2) Number of occurences of A+B / Number of occurences of B >= 50%
    # i.e., the conditional probability P(AB | A or B) > 50%

    # We have two dictionaries that keep track of our Pfam clusters.

    # ClusteringGroups assigns a number to each clustering of Pfams. This is an easy way to keep Pfams clumped up correctly.
    # PfamToClusteringID points each Pfam UID to its clustering group number. This saves us from having to search through
    # entries in the ClusteringGroups dictionary
    ClusteringGroups = {}
    PfamToClusteringID = {}

    for Pfam in PfamGroupingDictionary:
        TotalOccurences_Pfam = PfamGroupingDictionary[Pfam]['Total Alone'] + PfamGroupingDictionary[Pfam]['Total Combined']
        # We look at every linked pfam that shows up with the Pfam we're considering
        for LinkedPfam in PfamGroupingDictionary[Pfam]:
            # We ignore entries that are not Pfam IDs
            if LinkedPfam != 'Total Alone' and LinkedPfam != 'Total Combined':
                TotalOccurences_LinkedPfam = PfamGroupingDictionary[LinkedPfam]['Total Alone'] + PfamGroupingDictionary[LinkedPfam]['Total Combined']
                # P(AB | A)
                ProbABGivenA = PfamGroupingDictionary[Pfam][LinkedPfam]/TotalOccurences_Pfam*100

                # P(AB | B)
                ProbABGivenB = PfamGroupingDictionary[LinkedPfam][Pfam]/TotalOccurences_LinkedPfam*100
                # If we decide a link should be made, we add the pfams to the same clustering group
                if ProbABGivenA >= 50 or ProbABGivenB >= 50:

                    # If neither of the pfams we're considering are in any of the dictionaries yet, we add them
                    if Pfam not in PfamToClusteringID and LinkedPfam not in PfamToClusteringID:
                        if len(ClusteringGroups) != 0:
                            MaxPfamClusteringID = max(ClusteringGroups)
                        else:
                            MaxPfamClusteringID = 0
                        PfamClusteringID = MaxPfamClusteringID + 1
                        ClusteringGroups[PfamClusteringID] = [Pfam,LinkedPfam]
                        PfamToClusteringID[Pfam] = PfamClusteringID
                        PfamToClusteringID[LinkedPfam] = PfamClusteringID
                    # Otherwise we determine which one is and use the homology group number associated with it
                    # Note: we can do this because we only ever consider a pair once, so if a pfam ID is already
                    # in the dictionary, its pair won't be (this is because singlets are not considered)
                    else:
                        if Pfam in PfamToClusteringID:
                            PfamClusteringID = PfamToClusteringID[Pfam]
                            PfamToClusteringID[LinkedPfam] = PfamClusteringID
                            ClusteringGroups[PfamClusteringID].append(LinkedPfam)
                        else:
                            PfamClusteringID = PfamToClusteringID[LinkedPfam]
                            PfamToClusteringID[Pfam] = PfamClusteringID
                            ClusteringGroups[PfamClusteringID].append(Pfam)

                # If no link is made, the two pfams are given separate entries in the clustering group dictionary
                # If they already exist in the dictionary, nothing is done, otherwise the next clustering group number 
                # is generated by adding 1 to the highest existing 
                else:
                    if len(ClusteringGroups) == 0:
                        MaxPfamClusteringID = 0
                    else:
                        MaxPfamClusteringID = max(ClusteringGroups)
                    PfamClusteringID =  MaxPfamClusteringID + 1
                    if LinkedPfam not in PfamToClusteringID:
                        PfamToClusteringID[LinkedPfam] = PfamClusteringID
                        ClusteringGroups[PfamClusteringID] = [LinkedPfam]
                        PfamClusteringID += 1
                    if Pfam not in PfamToClusteringID:
                        PfamToClusteringID[Pfam] = PfamClusteringID
                        ClusteringGroups[PfamClusteringID] = [Pfam]
                        PfamClusteringID += 1
                # So that we don't look at the same pair twice (speeding up our algorithm), we delete the pfam that
                # we're searching in the outer for loop from its connected entry
                PfamGroupingDictionary[LinkedPfam].pop(Pfam)

    # The results are then dumped to flat files for use later
    json.dump(ClusteringGroups,open(ClusteringID_To_PfamCluster_Dictionary,'w'))
    json.dump(PfamToClusteringID,open(Pfam_To_ClusteringID_Dictionary,'w'))

########################################################################
#                        Program Executes Below                        #                     
########################################################################

# The user's table information is loaded into a dicitonary so we can loop through the tables
# If another database needs to be added, add a new column above and insert those entries
# into the dictionary below using the same format provided
DataTables = {'Ensembl' : {'ProteinTable' : EnsemblProteinTable,
                           'ProteinMetricsTable' : EnsemblProteinMetricsTable,
                           'UIDColumn': Ensembl_UIDColumn,
                           'PfamColumn': Ensembl_PfamColumn,
                           'HomGpColumn': Ensembl_HomGpColumn,
                           'ProteinMetricsUIDColumn': Ensembl_ProteinMetricsUIDColumn},
              
              'NCBI' : {'ProteinTable' : NCBIProteinTable,
                        'ProteinMetricsTable' : NCBIProteinMetricsTable,
                        'UIDColumn': NCBI_UIDColumn,
                        'PfamColumn': NCBI_PfamColumn,
                        'HomGpColumn': NCBI_HomGpColumn,
                        'ProteinMetricsUIDColumn': NCBI_ProteinMetricsUIDColumn}}


startTime = datetime.datetime.now()
if Verbose == True:
    print('\n\nProgram Executing\nCurrent Time: %s'%startTime)
    sys.stdout.flush()
    
# A connection is established to mysql 
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

if Verbose == True:
    print('Assembling Pfam Ages Dictionary\n')
    intermediateTime = datetime.datetime.now()
    sys.stdout.flush()
    
# We generate a dictionary where each pfam UID points to its age
PfamAgesDictionary = {}
mycursor.execute("SELECT " +','.join([PfamUIDColumn,AgeColumn]) + " FROM " + PfamDataTable )
PfamAgeResults = mycursor.fetchall()
for ageResult in PfamAgeResults:
    PfamUID,Age = ageResult
    PfamAgesDictionary[PfamUID] = float(Age)


# If the pfam clustering algorithm isn't being run, then the program attempts to
# load results from a previous run. If the files aren't found, then the user is
# notified and the program exits
if GeneratePfamClusters == False:
    if Verbose == True:
        print('Loading Pfam Clustering Dictionaries')
        sys.stdout.flush()
        
    try:
        PfamHomologyGroupLocation = json.load(open(Pfam_To_ClusteringID_Dictionary,'r'))
        PfamClusteringGroups = json.load(open(ClusteringID_To_PfamCluster_Dictionary,'r'))
        PfamClusteringGroups = {int(i):j for i,j in PfamClusteringGroups.items()}
    except:
        print('\n\nNo Pfam Clustering Dictionaries Found\nCheck Filenames or Set GeneratePfamClusters to True.\n\n')
        sys.exit(0)

# Otherwise, the program generates the relevant files and loads them back into the program
else:
    GeneratePfamClusteringDictionaries(PfamAgesDictionary,DataTables,Verbose,Pfam_To_ClusteringID_Dictionary,ClusteringID_To_PfamCluster_Dictionary)
    PfamHomologyGroupLocation = json.load(open(Pfam_To_ClusteringID_Dictionary,'r'))
    PfamClusteringGroups = json.load(open(ClusteringID_To_PfamCluster_Dictionary ,'r'))
    PfamClusteringGroups = {int(i):j for i,j in PfamClusteringGroups.items()}
    if Verbose == True:
        print('Pfam Clusters Generated\nTime Taken: %s\nAssigning Full-Protein Homology Groups\n'%(datetime.datetime.now()-intermediateTime))
        intermediateTime = datetime.datetime.now()
        sys.stdout.flush()

# Now the program uses the clustering information to assign homology groups to the full proteins in the database
# It does this by using the following algorithm:
# For simplicity, let's consider two clustering IDs. Lets say a protein has any number of oldest Pfams associated with it,
# we then look at the clustering IDs associated with those pfam IDs and assign homology groups based on those
#
# If A and B are clustering IDs:
#   - if A = B, then this is trivial and we assign the homology group number associated with that clustering ID
#   - if A =/= B, then we generate a new homology group number associated with AB


for DataTableSource in DataTables:

    ProteinTable = DataTables[DataTableSource]['ProteinTable']
    ProteinMetricsTable = DataTables[DataTableSource]['ProteinMetricsTable']
    UIDColumn = DataTables[DataTableSource]['UIDColumn']
    PfamColumn = DataTables[DataTableSource]['PfamColumn']
    HomGpColumn = DataTables[DataTableSource]['HomGpColumn']
    ProteinMetricsUIDColumn = DataTables[DataTableSource]['ProteinMetricsUIDColumn']
    
    if Verbose == True:
        print('Extracting %s Data'%DataTableSource)
        sys.stdout.flush()
        
    mycursor.execute("SELECT "+','.join([UIDColumn,PfamColumn]) + " FROM " + ProteinTable)

    if Verbose == True:
        print('Data Extracted\nTime Taken: %s\nParsing %s Results\n'%(datetime.datetime.now()-startTime,DataTableSource))
        intermediateTime = datetime.datetime.now()
        sys.stdout.flush()

    ProteinResults = mycursor.fetchall()
    for proteinResult in ProteinResults:
        # This dictionary will be used to keep track of the oldest pfams for
        # a particular protein and so will be reset for each entry
        proteinSpecificPfams = {}

        UID, Pfams = proteinResult
        PfamList = Pfams.split(',')

        for Pfam in PfamList:
            if Pfam in PfamAgesDictionary:
                AgeOfPfam = PfamAgesDictionary[Pfam]
                if AgeOfPfam not in proteinSpecificPfams:
                    # We make sure we only have one copy for each pfam that shows up
                    # in the protein by making the collection of PfamUIDs a set
                    proteinSpecificPfams[AgeOfPfam] = set([Pfam])
                else:
                    proteinSpecificPfams[AgeOfPfam].add(Pfam)
        # Ignores proteins without valid pfams
        if len(proteinSpecificPfams) != 0:
            # We then select the oldest age that shows up in a protein
            MaxAge = max(proteinSpecificPfams)
            # And collect the pfams that all share that common, maximum age
            MaxAgePfams = proteinSpecificPfams[MaxAge]

            # This should be cleaned up in the submodule at some point, but pfams that show up as the
            # oldest and don't share their age with any other pfam don't wind up in the clustering dictionaries
            # We add them below
            for MaxAgePfam in MaxAgePfams:
                # First we check if the pfam exists in the dictionary
                try:
                    PfamHomologyGroupLocation[MaxAgePfam]
                # If it doesn't, we assign it the next highest clustering ID in the dictionary
                except KeyError:
                    ClusteringID = max(PfamClusteringGroups) + 1

                    PfamHomologyGroupLocation[MaxAgePfam] = ClusteringID
                    PfamClusteringGroups[ClusteringID] = MaxAgePfam
            # Next, since the only thing we're interested in is the clusteringIDs that show up in each protein, we
            # reduce the set we're considering from the max age pfams to just the set of clustering IDs
            homologyGroupIDs = set([PfamHomologyGroupLocation[i] for i in MaxAgePfams])
            # If only one clustering ID is associated with that protein, then the homology group number assignment is
            # trivial
            if len(homologyGroupIDs) == 1:
                AssignedHomologyGroupID = list(homologyGroupIDs)[0]
            # If there's more than one, then we assign a new homology group number
            else:
                # It's possible that that combination of clustering IDs has shown up before, so we first to check to see
                # whether there's a homology ID associated with that cluster
                clustering = ','.join([str(j) for j in sorted(list(set([PfamHomologyGroupLocation[i] for i in MaxAgePfams])))])
                # If it has never shown up before, then we add it to our dictionaries and assign the new homology group ID
                # to that cluster
                if clustering not in PfamHomologyGroupLocation:
                    AssignedHomologyGroupID = max(PfamClusteringGroups) + 1
                    PfamHomologyGroupLocation[clustering] = AssignedHomologyGroupID
                    PfamClusteringGroups[AssignedHomologyGroupID] = clustering
                # If that particular cluster has shown up before, then we assign the homology group ID that's associated with that
                # cluster to the protein
                else:
                    AssignedHomologyGroupID = PfamHomologyGroupLocation[clustering]
            # Once a homology group ID has been assigned to the protein, that number is uploaded into the protein metrics dictionary
            UpdateStatement = "UPDATE " + ProteinMetricsTable + " SET " + HomGpColumn + "=" + str(AssignedHomologyGroupID) + " WHERE " + ProteinMetricsUIDColumn + "=" + str(UID)
            mycursor.execute(UpdateStatement)
            cnx.commit()

cnx.close()

if Verbose == True:
    print('All Proteins Successfully Assigned Homology Groups\nTime Taken: %s\nCurrent Time: %s'%(datetime.datetime.now()-startTime, datetime.datetime.now()))
    sys.stdout.flush()
