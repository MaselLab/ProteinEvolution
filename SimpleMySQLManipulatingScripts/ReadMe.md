

# Simple MySQL Manipulating Scripts
--------------------------------------

The MySQL database PFAMphylostratigraphy contains fairly massive data tables. This makes updating, transferring, and backing up tables complicated and, in some cases, impossible using MySQL GUIs such as Navicat.

This repository contains scripts that are designed to perform these simple tasks. 

This ReadMe is set up to give an overview of all the scripts contained in this repository. All scripts are listed below in alphabetical order.

In terms of editing tables and figuring out the computational/memory costs of operations and the consequences of making table alterations, I find [this](https://dev.mysql.com/doc/refman/5.6/en/innodb-online-ddl-operations.html) to be a good guide. For small SQL tables, oftentimes performing operations such as adding a column or changing the data type is quick and easy. However, it may or not be apparent to the user that it requires the table to be rebuilt, so performing this operation on gigantic datasets can take quite a long time and/or be impossible!

**Note: if a script is written in Python, it is likely necessary to acquire the mysql.connector module to be able to run. This can be easily installed from the command line with conda using the command found [here](https://anaconda.org/anaconda/mysql-connector-python). To be able to use conda, the user will need to have [anaconda3](https://www.anaconda.com/download/) installed.**

### AddIndexToMySQLTable.py
--------------
This script's function is to add an index to a MySQL table for a preexisting column. This is a really good idea to do for large tables where the user wishes to extract/upload entries based on values not associated with the table's primary key. Examples of this include:

   - Downloading a single species from a large, multi-species database by specifying species UID or species name
   - Uploading entries based on a key from another database. For example, in the PFAMphylostratigraphy databases, there are protein data tables that have their own keys as well as listing the coding table keys. If the user were to perform an analysis on a coding sequence and try to update the protein table without indexing the coding table column in the protein database, the upload process would take a very, very long time. 

Indexes can be added in Navicat by choosing the "design table" option, but this can be overly-slow and cumbersome for particularly large databases. This script will allow the user to add an index while letting it run in the background.

The script takes the following user-provided input:

   1) Database   : The name of the MySQL database where the relevant table is stored
   2) User       : The username to access MySQL
   3) Host       : The IP address to access MySQL
   4) Password   : The password to access MySQL
   5) DataTable  : The table where the user wants to add an index
   6) ColumnName : The column where the user wants to add the index



### BackupMySQLDataTable.py
--------------
This script is written using python3 and is intended to be used with very large data tables that cannot be backed up using more conventional methods. For example, the data table on MySQL in the database PFAMphylostratigraphy called NCBIGenomes_Protein_Complete has ~seven-million entries and cannot be backed up using the copy command in Navicat or the DUPLICATE SQL command.

This script takes as input the following:

  - [ ] Database - This is the name of the database where the data table being backed up is located
  - [ ] User - The username to access MySQL
  - [ ] Host - The IP address to access MySQL
  - [ ] Password - User's password to access MySQL
  - [ ] TableToBeBackedUp - This is the table the user is backing up
  - [ ] NewTableName - This is the name of the table that will be created as the backup table
  
Once running, the script performs a check to make sure a table with the new table name doesn't already exist. This is to prevent the inadvertent overwriting of data.

The script backs up by pulling the entries from the source table one-by-one and inserting them into the new table. It commits these entries in batches of 1,000. 

Time: In general, this script uploads and commits entries at a rate of approximately 1,000 entries/second. For a table like NCBIGenomes_Protein_Complete, this amounts to a runtime of roughly 2.5 hours. As a result, the user may wish to run this script with nohup 
