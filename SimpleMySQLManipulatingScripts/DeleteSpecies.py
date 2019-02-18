import os, json, sys, csv, mysql.connector, datetime

'''
Author : Sara Willis
Date   : February 18, 2019

This is a very simple script that's used to delete a specific species from a mysql database using a species UID
'''

DataTable = 'TestUploadTable_Protein_UseMeForTests'
SpeciesUID = 294

# User-specific MySQL connection information
Database = 'PFAMphylostratigraphy'
User = ''
Host = ''
Password = ''

# Establishes a connection to MySQL
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Defines a deletion statement and tuple with the species UID to delete
DeleteSpeciesStatement = "DELETE FROM "+DataTable+" WHERE SpeciesUID = %s"
DeletionData = tuple((SpeciesUID,))

# Executes the deletion statement and commits the changes
mycursor.execute(DeleteSpeciesStatement,DeletionData)
cnx.commit()
cnx.close()
