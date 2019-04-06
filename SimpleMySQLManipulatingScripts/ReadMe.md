

# Simple MySQL Manipulating Scripts
--------------------------------------

The MySQL database PFAMphylostratigraphy contains fairly massive data tables. This makes updating, transferring, and backing up tables complicated and, in some cases, impossible using MySQL GUIs such as Navicat.

This repository contains scripts that are designed to perform these simple tasks. 

This ReadMe is set up to give an overview of all the scripts contained in this repository. All scripts are listed below in alphabetical order.

In terms of editing tables and figuring out the computational/memory costs of operations and the consequences of making table alterations, I find [this](https://dev.mysql.com/doc/refman/5.6/en/innodb-online-ddl-operations.html) to be a good guide. For small SQL tables, oftentimes performing operations such as adding a column or changing the data type is quick and easy. However, it may or not be apparent to the user that it requires the table to be rebuilt, so performing this operation on gigantic datasets can take quite a long time and/or be impossible!

If you need any guidance with respect to executing MySQL commands, [this is a good resource](https://www.w3schools.com/sql/sql_delete.asp).

**Note: if a script is written in Python, it is likely necessary to acquire the mysql.connector module to be able to run. This can be easily installed from the command line with conda using the command found [here](https://anaconda.org/anaconda/mysql-connector-python). To be able to use conda (and Python3 on Fusion), the user will need to have [anaconda3](https://www.anaconda.com/download/) installed.**

# Python Scripts

### AddIndexToMySQLTable.py

Author : Sara Willis
Date   : Monday January 7, 2019

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

### AddMySQLColumn

Author : Sara Willis
Date   : March 11, 2019

This script is designed to add a column to a MySQL table. When working with relatively small tables, this can be done using the GUI Navicat. As tables get larger this becomes a more cumbersome because when a column is added, it requires the table be recreated which can take a long time. This script allows Navicat to remain in an unfrozen state while a column is added to a table in the background.

### DeleteSpecies

Author : Sara Willis
Date   : February 18, 2019

This is a very simple script that's used to delete a specific species from a mysql database using a species UID

### MakeDataTableBackup.py

Author: Sara Willis
Date  : Monday January 7, 2019

This script is used to make backups of large MySQL data tables that cannot be backed up by more conventional methods
This script pulls entries from the table that needs to be backed up and inserts them into a target table, committing them in batches of 1000.

-----------------DEPENDENCIES----------------

The user will need the [python module mysql.connector](https://anaconda.org/anaconda/mysql-connector-python) for this script to run 


-------------------RUNTIME-------------------

This process can take some time depending on the size of the data table that is being backed up. The user may wish to run this process with nohup. For the table NCBIGenomes_Protein_Complete, this script can take about 2.5 hours to run through.
