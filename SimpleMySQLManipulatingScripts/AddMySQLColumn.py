import os, sys, mysql.connector, datetime

'''

Author : Sara Willis
Date   : March 11, 2019


This script is designed to add a column to a MySQL table. When working with relatively small tables, this can be done using the GUI Navicat. As tables get larger this becomes a more cumbersome because when a column is added, it requires the table be recreated which can take a long time. This script allows Navicat to remain in an unfrozen state while a column is added to a table in the background. 

'''

####################################################################################
#                                  User Options                                    #
####################################################################################


Database = ''                  
User = ''                                     
Host = ''                                      
Password = ''

DataTable = 'Genomes_Multicellular_ProteinMetrics'         
ColumnToBeAdded = 'DeltaTango_DensityOfAggProneRegions'  
PrecedingColumn = 'DensityOfAAsInAPRs'     # If this is left as blank (''), the column will automatically be added as the last
                            # If this argument is set to 'FIRST', the column will be added as the first
                            # If an existing column name is given, then the new one will be added immediately following it
                      
ColumnType = 'DECIMAL'      # Options: INT, FLOAT, TEXT, VARCHAR, DECIMAL
ColumnLength = '5'          # If this is left blank, then the column defaults to the largest allowable value

NumberOfDecimals = '4'      # Only used when ColumnType = DECIMAL


####################################################################################
#                             Program Executes Below                               #
####################################################################################


start_time = datetime.datetime.now()

# Establish a connection to the mysql database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Sets up the column length. It's unneccessary to add anything if you want the max length available for the given column type
# otherwise, if only one length is specified, such as for INT, the length should be included with the variable as VARIABLE(length)
# If the variable type is DECIMAL, it requires two values, the length of the decimal and the number of decimals: DECIMAL(length, number of decimals)
# For example:
#    xx.xxx = DECIMAL(5,3)
#    x.xxxx = DECIMAL(5,4)
if ColumnLength != '':
    if NumberOfDecimals != '':
        LengthVariable = '(%s,%s)'%(ColumnLength,NumberOfDecimals)
    else:
        LengthVariable = '(%s)'%ColumnLength
else:
    LengthVariable = ''

# Here we specify where in the table the included column goes
if PrecedingColumn == '':
    LocationVariable = ''
elif PrecedingColumn == 'FIRST':
    LocationVariable = 'FIRST'
else:
    LocationVariable = 'AFTER %s'%PrecedingColumn

# And the column is added
AddColumnStatement = "ALTER TABLE "+DataTable+" ADD COLUMN "+ColumnToBeAdded+" " +ColumnType+LengthVariable + ' ' +LocationVariable
mycursor.execute(AddColumnStatement)
cnx.commit()
cnx.close()
print('\n\nColumn Successfully Added\nTime Taken: %s\n\n' %(datetime.datetime.now()-start_time))


