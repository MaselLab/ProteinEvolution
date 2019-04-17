import mysql.connector

'''
Author : Sara Willis
Date   : April 17, 2019


The purpose of this script is to calculate the size of a MySQL datatable in MB for the user. The user should enter their MySQL connection information and the name of the table they want to know the size of under User-Specific Data. The program will then print the size of the table to the terminal 
'''


##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################

# MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# Name of the table to determine the size of 
TableName = 'NCBIGenomes_Protein_Complete'

##########################################################################################
#                                Program Executes Below                                  #
##########################################################################################

# Program establishes a connection to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We define the statement that will pull the size of the table in MB
MySQLStatement = "SELECT table_name AS 'Table', round(((data_length + index_length)/1024/1024),2) 'Size in MB' FROM information_schema.TABLES WHERE table_schema = '"+Database+"' AND table_name = '"+TableName+"'"

mycursor.execute(MySQLStatement)

# Prints the size of the table to the terminal for the user
print('\nSize of Table "'+TableName+'" : %s MB\n'%mycursor.fetchall()[0][1])

    
