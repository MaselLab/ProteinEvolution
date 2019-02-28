import os, json, sys, csv, mysql.connector, datetime, copy
import numpy as np

'''
Author : Sara Willis
Date   : February 11, 2019
--------------------------


In each protein data table in the database PFAMphylostratigraphy, there is one row for for each protein and its associated Pfams. These tables have then been deconstructed to make the DomainMetrics data tables where there is one entry for each protein/pfam combo. What this script does is it gathers the UIDs in the DomainMetrics data tables, associates all the UIDs that belong to a particular protein, and uploads those back into the protein data table. This allows for unique keys that allow the user to go back and forth between data tables without ambiguity.

============================================================
Example of Script's Function:

If a protein is stored as the following in a protein data table:

------------------------------------------------------------
ProteinUID  |           PfamUIDs      
____________|_______________________________________________
    1       |   PF00001,PF00002,PF00003
------------------------------------------------------------

Then it has three pfams associated with it that are stored with the following format in the DomainMetrics data table:

------------------------------------------------------------
PfamTableUID   |    ProteinUID    |   PfamUID
_______________|__________________|_________________________
     1         |        1         |   PF00001
     2         |        1         |   PF00002
     3         |        1         |   PF00003
------------------------------------------------------------

What this script accomplishes is it inserts the PfamTableUIDs associated with each protein back into the protein data table so it's easy to move back and forth between the tables. The final format of the protein table should be the following:

------------------------------------------------------------
ProteinUID  |           PfamUIDs            | PfamTableUIDs
____________|_______________________________|_______________
    1       |   PF00001,PF00002,PF00003     |     1,2,3
------------------------------------------------------------
'''

# The user should enter their login information to access MySQL here
Database = 'PFAMphylostratigraphy'
User = ''
Host = ''
Password = ''

DomainMetricsTable = 'EnsemblGenomes_DomainMetrics_Complete'
ProteinTable = 'EnsemblGenomes_Protein_Complete'

# The script then connects to the database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We select the MaxUID so that we can iterate through the table one row at a time
ExtractionStatement = "SELECT MAX(UID) FROM "+DomainMetricsTable
mycursor.execute(ExtractionStatement)
MaxUID = mycursor.fetchone()[0]

'''
The general principal is the following: A for loop was used to the generate the domain metrics data tables, iterating through each protein in the source data table, and then iterating through each PfamUID associated with the protein, making a row for each PfamUID/protin combo. This algorithm reconstructs this logic. It goes through each entry in the domain metrics data table and start gathering UIDs. It collects for each protein UID in the table until the protein UID changes. Once the protein UID changes, we've found all Pfam UIDs associated with that protein in the domains data table. We take all the UIDs we've found, concatenate them into a comma-delimited string and upload them into the protein data table. We then start the process again. We do this until we've exhausted the entire table
'''

# Start tells us we're beginning the search. This is so we can get our variables set up to use throughout
# the duration of the script
Start = True

# We iterate through each UID in the protein table and add one to the max UID found because the range(1,n) function
# goes up to n-1 by design
for i in range(1,MaxUID+1):
    SelectUIDDAndDomain = "Select UID,ProteinTableUID FROM " + DomainMetricsTable + " WHERE UID = %s"%i
    mycursor.execute(SelectUIDAndDomain)
    result = mycursor.fetchone()
    if result != None:
        UID,ProteinTableUID = result
                
        if Start == True:
            # We keep track of the last protein UID that we found. This way we can determine when we've moved
            # on to a new protein 
            lastProteinUID = copy.deepcopy(ProteinTableUID)
            Protein_Pfam_Correspondance = [UID]
            Start = False
        else:
            # If this is just a continuation of a previous protein, we add the pfam UID to our list and move on
            if ProteinTableUID == lastProteinUID:
                Protein_Pfam_Correspondance.append(UID)
            else:
                # If, however, our protein has changed and we've found everything for the previous protein
                # we need to upload our results for the last protein we searched. We join all the UIDs for that
                # protein into a comma-delimited string and use them to update the protein table
                PfamUIDs = ','.join([str(j) for j in Protein_Pfam_Correspondance])
                UpdateStatement = "UPDATE "+ProteinTable+" SET PfamMetricsTableUID=%s WHERE UID = %s"
                UpdateData = (PfamUIDs, lastProteinUID)
                mycursor.execute(UpdateStatement,UpdateData)
                cnx.commit()
                # Once the last protein has been dealt with, we start a new list for the next protein we've found
                Protein_Pfam_Correspondance = [UID]
            lastProteinUID = copy.deepcopy(ProteinTableUID)

# Once we've run out of entries from our table, we upload the results for the last protein in our data table, and
# close the connection.
UpdateStatement = "UPDATE "+ProteinTable+" SET PfamMetricsTableUID=%s WHERE UID = %s"
UpdateData = (PfamUIDs, lastProteinUID)
mycursor.execute(UpdateStatement,UpdateData)
cnx.commit()
