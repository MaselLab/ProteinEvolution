import os, sys, json, csv, mysql.connector, datetime



'''

Created :  Monday January 7, 2019

Author  :  Sara Willis



This script is used to make backups of large MySQL data tables that cannot be backed up by more conventional methods

This script pulls entries from the table that needs to be backed up and inserts them into a target table, committing them in batches of 1000.

-----------------DEPENDENCIES----------------

The user will need the python module mysql.connector for this script to run 

https://anaconda.org/anaconda/mysql-connector-python



-------------------RUNTIME-------------------
This process can take some time depending on the size of the data table that is being backed up. The user may wish to run this process with nohup. For the table NCBIGenomes_Protein_Complete, this script can take about 2.5 hours to run through. 
'''



###############################################################
######## Enter the user-specific information below ############
###############################################################


# The user needs to enter their information to access the MySQL database
Database = ''  # Database where data table is located
User = ''      # Username to access MySQL 
Host = ''      # IP address to access MySQL
Password = ''  # Password to access MySQL

# The name of the table that needs to be backed up is specified here
TableToBeBackedUp = 'TestUploadTable_Protein_UseMeForTests'

# And the name of the new table where the data will be copied also needs to be specified
NewTableName = 'TestUploadTable_Protein_UserMeForTests_copy'

# Verbose, when set to True, notifies the user of the programs progress and how long the operation has taken
Verbose = True

###############################################################






###################################################################################################
#################################### PROGRAM EXECUTES BELOW #######################################
###################################################################################################

# Starts timing the program:
start_time = datetime.datetime.now()

# Commits 1000 entries per batch 
batchIntervalSize = 1000

# Connects to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)

# Defines a cursor so the program can interact with the database
mycursor = cnx.cursor(buffered = True)

# We want the column names for the data in the table we want to back up. We do this by constructing
# a list containing all the column names. Later we'll modify this to a string to be used in an
# insertion statement 
mycursor.execute("SHOW COLUMNS FROM %s"%TableToBeBackedUp)
Columns = mycursor.fetchall()
ColumnList = []
for column in Columns:
    ColumnName = column[0]
    ColumnType = column[1]
    ColumnList.append(ColumnName)

# First, to make sure that the user hasn't chosen to name their table after an already-existing table
# we try to select an entry from the table name specified. This is to prevent involuntarily overwriting
# a preexisting table with potentially catastrophic consequences! Be careful with your data!
try:
    
    # If we successfully extract an entry, then we know the table already exists. If this is the case,
    # the user is notified and the program exits.
    mycursor.execute("SELECT 1 FROM %s LIMIT 1" %NewTableName)
    print('\n\nThe new table name you have specified already exists!\nPlease choose a valid table name and try again\n\n')
    
    
except:
    
    # If the program was unable to extract anything from the table the user has designated as the backup,
    # then the table does not exist and the backup proceeds. We start by duplicating the structure of
    # the table we want to back up.
    if Verbose == True:
        print('\n\nCreating Backup Table\nCurrent Time: %s\n\n'%(datetime.datetime.now()))
        sys.stdout.flush()
    mycursor.execute("CREATE TABLE %s LIKE %s"%(NewTableName,TableToBeBackedUp))
    cnx.commit()

    # Pulls the max UID in the table that's being backed up.
    # This is used so we can run through the entire table, pulling
    # only one entry at a time. This prevents us from having to pull all entries at once which
    # can take some time, depending on the table size.
    
    GetMaxUID = 'SELECT MAX(UID) FROM %s'%TableToBeBackedUp
    mycursor.execute(GetMaxUID)
    MaxUID = mycursor.fetchone()[0]

    # Runs through all UIDs in the table starting at 1 and stopping at the MaxUID
    # (1 is added since range(1,n) is of the form (1,...,n-1))
    
    for n in range(1,MaxUID+1):
        
        # Each row is extracted using n to pull entry with UID = n
        extractProteinSequence = 'SELECT * FROM %s WHERE UID = %s'%(TableToBeBackedUp,n)
        mycursor.execute(extractProteinSequence)
        
        # Since only one entry was pulled, we use fetchone so no for loop needs to be run
        # to retrieve entries. This consolidates the results into a single tuple as opposed
        # to a list of tuples
        result = mycursor.fetchone()
        
        # If an entry with UID = n doesn't exist, then nothing is done
        # This can happen if entries were removed from the original table.
        if result == None:
            pass
        # If the entry does exist, then it is inserted into the new table
        else:
            # An insertion statement is defined. We do this by constructing multiple parts:
            #   1) The insertion command INSERT INTO
            #   2) The new table name where we're copying our data NewTableName
            #   3) The list of columns reformatted as a string so SQL can read it
            #   4) VALUES -- the values that will be inserted follow this statement
            #   5) And a string of the form (%s,...,%s) so we can insert our custom values into the table
            sqlUploadProteinDataStatement = 'INSERT INTO '+ NewTableName +' '+(str(tuple(ColumnList)).replace("'",'')) +' VALUES ' + '('+','.join(['%s' for i in ColumnList]) + ')'

            # And the row is inserted using the execute command
            mycursor.execute(sqlUploadProteinDataStatement,result)
            
            # If n exceeds the maxUID value in the table being copied,
            # the upload stops, the entries are committed and the
            # program shuts down
            if n > MaxUID:
                cnx.commit()
                cnx.close()
                sys.exit(0)
            # Each interval of 1000 entries are committed to the table
            if n>batchIntervalSize:
                batchIntervalSize += 1000
                cnx.commit()

    # Prints the time it took for the backup to be made
    if Verbose == True:
        print('Backup Complete!\nDuration: %s\n\n'%(datetime.datetime.now()-start_time))
        sys.stdout.flush()

    # Once the table is completely copied, the last of the entries are committed
    cnx.commit()
    cnx.close()


