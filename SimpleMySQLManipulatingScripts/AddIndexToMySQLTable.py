import json, sys, os, csv, mysql.connector, datetime


'''
The following script was created on Monday January 7, 2019

Author: Sara Willis


This script's function is to add an index to a MySQL table for a preexisting column. This is a really good idea to do for large tables where the user wishes to extract/upload entries based on values not associated with the table's primary key. Examples of this include:

   1) Downloading a single species from a large, multi-species database by specifying species UID or species name
   2) Uploading entries based on a key from another database. For example, in the PFAMphylostratigraphy databases, there are protein data tables that have their own keys as well as listing the coding table keys. If the user were to perform an analysis on a coding sequence and try to update the protein table without indexing the coding table column in the protein database, the upload process would take a very, very long time. 

Indexes can be added in Navicat by choosing the "design table" option, but this can be overly-slow and cumbersome for particularly large databases. This script will allow the user to add an index while letting it run in the background.

The script takes the following user-provided input:

   Database   -- The name of the MySQL database where the relevant table is stored
   User       -- The username to access MySQL
   Host       -- The IP address of the MySQL database
   Password   -- The password to access MySQL
   DataTable  -- The table the user wants to add an index to
   ColumnName -- The column the user wants to add the index to



DEPENDENCIES:
-------------

mysql.connector: https://anaconda.org/anaconda/mysql-connector-python

'''


#----------    User Information    ----------#

Database = ''    # The name of the MySQL database where the relevant tables are stored
User = ''        # The username to access MySQL
Host = ''        # The IP address of the MySQL database
Password = ''    # The password to access MySQL
DataTable = ''   # The table where the user wants to add an index
ColumnName = ''  # The column where the index is being added. Verify that this is a string and not an integer. If it's an integer, it's possible there will be some unforseen consequences.


#---------- Program Executes Below ----------#

# The program starts by setting up a connection with the database
cnx = mysql.connector.connect(user = User,
                            password = Password,
                            host = Host,
                            database = Database)

# It then creates a cursor so it's able to interact with the database
mycursor = cnx.cursor(buffered = True)

# A statement is defined to add an index to the data table
CreateIndexStatement = "ALTER TABLE %s ADD INDEX %s (%s)" %(DataTable,ColumnName,ColumnName)

# We then attempt to execute the statement
try:
    mycursor.execute(CreateIndexStatement)
    cnx.commit()
    cnx.close()
    
# If the execution of the statement doesn't work, the field already is indexed
except:
    print('\n\nIndex already exists for listed column!\n\n')
    cnx.close()
    

