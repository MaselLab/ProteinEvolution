

# SimpleMySQLManipulatingScripts
--------------------------------------


The MySQL database PFAMphylostratigraphy contains fairly massive data tables. This makes updating, transferring, and backing up tables complicated and, in some cases, impossible using MySQL GUIs such as Navicat.

This repository contains scripts that are designed to perform these simple tasks.


## BackupTable.py
--------------
This script is written using python3 and is intended to be used with very large data tables that cannot be backed up using more conventional methods. For example, the data table on MySQL in the database PFAMphylostratigraphy called NCBIGenomes_Protein_Complete has ~seven-million entries and cannot be backed up using the copy command in Navicat or the DUPLICATE SQL command.

This script takes as input the following:

  - [ ] Database - This is the name of the database where the data table being backed up is located
  - [ ] User - The username to access the mysql database
  - [ ] Host - This is the IP address to access mysql
  - [ ] Password - The user's password to access the mysql database
  - [ ] TableToBeBackedUp - This is table the user is backing up
  - [ ] NewTableName - This is the name of the table that will be created
  
Once running, the script performs a check to make sure a table with the new table name doesn't already exist. This is to prevent the inadvertent overwriting of data. If a table with that name is found, the script notfies the user and shuts down. Otherwise, it proceeds.

The script backs up by pulling the entries from the source table one-by-one and inserting them into the new table. It commits these entries in batches of 1,000. 

Time: In general, this script uploads and commits entries at a rate of approximately 1,000 entries/second. For a table like NCBIGenomes_Protein_Complete, this amounts to a runtime of roughly 2.5 hours. As a result, the user may wish to run this script with nohup 
