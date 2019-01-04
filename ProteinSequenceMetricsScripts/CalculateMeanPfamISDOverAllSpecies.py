'''

Created Friday January 4, 2019

Author: Sara Willis

The purpose of this script is to calculate the means and variances of all Pfam ISD values for all occurrences over all species in the PFAMphylostratigraphy database. The values are calculated by extracting from the metrics datatables:

    1) EnsembleGenomes_Metrics_Complete
    2) NCBIGenomes_Metrics_Complete

Once these values are calculated, they are uploaded into the data table PfamUIDsTable_EnsemblAndNCBI

'''



import os, sys, json,csv, datetime,mysql.connector
import numpy as np


# The time taken to execute this script will be returned to the user 
start_time = datetime.datetime.now()

# Following each major block, the user is notified of the time taken and which process has been and is about to be run
print('Starting Analysis\nCurrent Time: %s'%start_time)

# First, we connect to the database with the user's information
def LoadMySQLConnectionInformation():
    Database = ''
    User = ''
    Host = ''
    Password = ''
    return Database, User, Host, Password
Database, User, Host, Password = LoadMySQLConnectionInformation()
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
# We define a cursor so we can interact with the tables
mycursor = cnx.cursor(buffered=True)

# We start by extracting the metrics from the Ensembl table
ExtractionStatement = "SELECT ProteinTableUID,PfamUID,MeanISD_FullProtein,MeanISD_PfamOnly,SpeciesUID FROM EnsemblGenomes_Metrics_Complete"
print('Extracting Ensembl Metrics\n')
mycursor.execute(ExtractionStatement)
results = mycursor.fetchall()
print('Ensembl Metrics Extracted\nTime Taken: %s\n\nCreating Dictionary\n'%(datetime.datetime.now()-start_time))
currentTime = datetime.datetime.now()

# We then define two dictionaries where we will store our data:
#    1) PfamISDBySpecies -- This will allow us to differentiate the mean ISD values for each Pfam by species. This can be dumped to a dictionary
#    2) PfamISDFullDataset -- This is where each Pfam value is stored. The ISD values for each occurrence of the Pfam will be stored in a list
PfamISDBySpecies = {}
PfamISDFullDataset = {}

# The results from the MySQL table extraction are then parsed
for result in results:
    ProteinTableUID, PfamUID, MeanISD_FullProtein,MeanISD_PfamOnly,SpeciesUID = result
    # Here we ensure that the ISD values are floats and not strings which would cause problems later
    MeanISD_FullProtein, MeanISD_PfamOnly = float(MeanISD_FullProtein), float(MeanISD_PfamOnly)

    # If the species UID is not in the species dictionary, we add it. It points to a subdictionary that contains each Pfam that shows up in that species and the ISD values associated with that Pfam stored as a list
    if SpeciesUID not in PfamISDBySpecies:
        PfamISDBySpecies[SpeciesUID] = {PfamUID :[MeanISD_PfamOnly]}
    # If the species is already in the dicionary, we check to see whether the current Pfam is already in that subdictionary
    else:
        # If it's not, it's added. The PfamUID points to a list containing the mean ISD value for that entry
        if PfamUID not in PfamISDBySpecies[SpeciesUID]:
            PfamISDBySpecies[SpeciesUID][PfamUID] = [MeanISD_PfamOnly]
        # If the Pfam is already in the subdictionary, we add the current ISD value to it's list
        else:
            PfamISDBySpecies[SpeciesUID][PfamUID].append(MeanISD_PfamOnly)

    # The dictionary that is used to update the MySQL datatable is simpler. It only has Pfam UIDs as keys pointing to the list of the ISD values associated with it

    # First we check to see whether the PfamUID is already in the dictionary. If it isn't it's added
    if PfamUID not in PfamISDFullDataset:
        PfamISDFullDataset[PfamUID] = [MeanISD_PfamOnly]
    # Otherwise, the current ISD value is added to its list
    else:
        PfamISDFullDataset[PfamUID].append(MeanISD_PfamOnly)


# The user is notified of the progress
print('Dictionary Created\nTime Taken: %s\n\nExtracting NCBI Metrics\n'%( datetime.datetime.now()-currentTime))
# And the intermediate time is updated
currentTime = datetime.datetime.now()

# Next, we extract the data from the NCBI database
ExtractionStatement = "SELECT ProteinTableUID,PfamUID,MeanISD_FullProtein,MeanISD_PfamOnly,SpeciesUID FROM NCBIGenomes_Metrics_Complete"
mycursor.execute(ExtractionStatement)
results = mycursor.fetchall()
print('NCBI Metrics Extracted\nTime Taken: %s\n\nUpdating Dictionary\n'%(datetime.datetime.now()-currentTime))
currentTime =  datetime.datetime.now()

# We then go through the exact same process as we did for the Ensembl database 
for result in results:
    ProteinTableUID, PfamUID, MeanISD_FullProtein,MeanISD_PfamOnly,SpeciesUID = result
    MeanISD_FullProtein, MeanISD_PfamOnly = float(MeanISD_FullProtein), float(MeanISD_PfamOnly)
    if SpeciesUID not in PfamISDBySpecies:
        PfamISDBySpecies[SpeciesUID] = {PfamUID :[MeanISD_PfamOnly]}
    else:
        if PfamUID not in PfamISDBySpecies[SpeciesUID]:
            PfamISDBySpecies[SpeciesUID][PfamUID] = [MeanISD_PfamOnly]
        else:
            PfamISDBySpecies[SpeciesUID][PfamUID].append(MeanISD_PfamOnly)
    if PfamUID not in PfamISDFullDataset:
        PfamISDFullDataset[PfamUID] = [MeanISD_PfamOnly]
    else:
        PfamISDFullDataset[PfamUID].append(MeanISD_PfamOnly)
        
#json.dump(open('PfamISDFullDataset_NCBIPlusEnsembl.dictionary','w'))

# The user is notified of the progress
print('Dictionary Updated\nTime Taken: %s\n\nUpdating Pfam Table\n'%(datetime.datetime.now()-currentTime))

# And the intermediate time is updated
currentTime = datetime.datetime.now()

# We then extract the Pfams and their UIDs from the Pfam Table
extractPfamStatement = "SELECT UID,PfamUID FROM PfamUIDsTable_EnsemblAndNCBI"
mycursor.execute(extractPfamStatement)
results = mycursor.fetchall()

# We then check each Pfam that shows up in that table
for result in results:
    UID, PfamUID = result
    # We check to make sure that the Pfam UID is in the dictionary we've created. If we just assume it does, then if we try to process a subset of our full dataset, we'll wind up quitting with an error
    if PfamUID in PfamISDFullDataset:
        # We then calculate the mean ISD of each Pfam and its variance
        PfamMeanISD = np.mean(PfamISDFullDataset[PfamUID])
        PfamVariance = np.var(PfamISDFullDataset[PfamUID])
        # We then update the MySQL table
        PfamTableUpdateStatement = "Update PfamUIDsTable_EnsemblAndNCBI SET ISDOverAllSpecies_PfamOnly_Mean=%s,ISDOverAllSpecies_PfamOnly_Variance=%s WHERE UID = %s"%(PfamMeanISD,PfamVariance,UID)
        mycursor.execute(PfamTableUpdateStatement)
        # And commit the changes
        cnx.commit()

# The user is notified that the script executed successfully and how long it took.        
print('Pfam Table Updated\nTime Taken: %s\n\nScript Complete!\nTotal Time Taken: %s'%((datetime.datetime.now()-currentTime),(datetime.datetime.now()-start_time)))
    
