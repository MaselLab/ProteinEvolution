import os, sys, mysql.connector, datetime

'''
Author : Sara Willis
Date   : May 8, 2019

This script is designed to drop unwanted columns from a MySQL table. The user needs to define the column(s) that need to be deleted and the table name. If more than one column needs to go, then the user should input the columns as comma-delimited strings:


 ColumnsToBeDeleted = 'Column1' [, 'Column2'[ , 'Column3' [, ...]]]

 Where options in the square brackets are optional

 for example:

   ColumnsToBeDeleted = 'Column1'

for one column, or:

   ColumnsToBeDeleted = 'Column1','Column2',...

for multiple

'''

####################################################################################
#                                  User Options                                    #
####################################################################################


Database = ''
User = ''
Host = ''
Password = ''

DataTable = 'Genomes_Protists_ProteinMetrics'

# Add the names of all columns to be deleted from specified table as comma-delimited strings:
ColumnsToBeDeleted = 'HomologyGroupID'


####################################################################################
#                             Program Executes Below                               #
####################################################################################

# The program notifies the user of the progress of the script
start_time = datetime.datetime.now()

# If only one column needs to be deleted, the variable ColumnsToBeAdded will be a string
# and needs to be made into a list so we can use the .join function. If multiple columns
# are specified, the variable is a tuple
if type(ColumnsToBeDeleted) == str:
    ColumnsToBeDeleted = [ColumnsToBeDeleted]

print('\n\nProgram Executing\nAltered Table: %s\nDropped Columns: '%DataTable+','.join(ColumnsToBeDeleted) + '\nCurrent Time: %s'%start_time)


# Establish a connection to the mysql database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)


# We define the relevant DROP statement to get rid of the user-defined columns from the relevant table   
DropColumnsStatement = "ALTER TABLE " + DataTable +','.join([' DROP COLUMN '+i for i in ColumnsToBeDeleted])

mycursor.execute(DropColumnsStatement)
cnx.commit()

print('Column(s) Successfully Dropped From Table\nTime Taken: %s\n\n'%(datetime.datetime.now()-start_time))

cnx.close()
