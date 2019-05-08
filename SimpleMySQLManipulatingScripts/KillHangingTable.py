import os, json, sys, csv, mysql.connector, datetime, copy, matplotlib, requests

'''
Author : Sara Willis
Date   : Wednesday May 8, 2019


This script is designed to kill MySQL processes that are "hanging". This can happen if a user attempts to alter a table manually (i.e. using NaviCat or some other GUI) and the program freezes up. This can result in the user force-quitting the application, opening the application again, and attempting to open the table only to find that nothing loads.

To solve this problem, you can execute a command in MySQL that shows all the processes that are currently running. The problematic processes will show as: 
      

           "Waiting for table metadata lock"

This script will collect all the processes that are currently running on MySQL and will print them out on the command line for the user. It will then prompt them for the ID of the process they want to terminate. If they select a valid ID, the script will terminate it.

NOTE: This script will only allow the user to see the process they are currently running and not other users logged in remotely to the machine hosting MySQL so if there is another user locking up a table you will need to talk them into using this script. 
'''
########################################################################
#                        User-Specific Information                     #                     
########################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''


########################################################################
#                        Program Executes Below                        #                     
########################################################################

# Establishes a connection to MySQL
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# The Program will delete as many processes as the user wants
DeleteProcesses = 'y'

while DeleteProcesses == 'y':

    # All running processes are extracted and shown to the user
    mycursor.execute("SHOW PROCESSLIST")

    processes = mycursor.fetchall()

    # The IDs are kept track of so the input UID from the user can be cross-referenced with existing process UIDs
    ProcessIDs = []

    # The user can tinker with these margins if the columns are too narrow/wide. The margins are the number of spaces
    # wide each column is
    IDMargin = 10
    UserMargin = 20
    HostMargin = 20
    DatabaseMargin = 25
    CommandMargin = 20
    TimeMargin = 10
    StateMargin = 20
    InfoMargin = 30

    # The program prints a handy table that makes it easy to visualize the processes and their information
    print(('\n\n{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format('='*IDMargin,'='*UserMargin,'='*HostMargin,'='*DatabaseMargin,'='*CommandMargin,'='*TimeMargin,'='*StateMargin,'='*InfoMargin))
    print(('{:^%s}{:^%s}{:^%s}{:^%s}{:^%s}{:^%s}{:^%s}{:^%s}'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format('ID','USER','HOST','DATABASE','COMMAND','TIME','STATE','INFO'))
    print(('{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format('-'*IDMargin,'-'*UserMargin,'-'*HostMargin,'-'*DatabaseMargin,'-'*CommandMargin,'-'*TimeMargin,'-'*StateMargin,'-'*InfoMargin))

    # Each process is printed to the table
    for process in processes:
        ID,Usr,ProcessHost,DB,Command,Time,State,Info,extra = process
        if Info == None:
            Info = ''

        # And it's ID is recorded
        ProcessIDs.append(ID)

        # Most of the mess is the formatting
        print(('{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format('|','|','|','|','|','|','|','|','|'))
        print(('{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format("|   "+str(ID),"|  "+str(Usr),"|  "+str(ProcessHost),"|  "+str(DB),"|  "+str(Command),"|  "+str(Time),"|  "+str(State),"|  "+str(Info)))

    # Ends the column with a line footer
    print(('{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}{:<%s}\n\n'%(IDMargin,UserMargin,HostMargin,DatabaseMargin,CommandMargin,TimeMargin,StateMargin,InfoMargin)).format('='*IDMargin,'='*UserMargin,'='*HostMargin,'='*DatabaseMargin,'='*CommandMargin,'='*TimeMargin,'='*StateMargin,'='*InfoMargin))

    # The program asks the user if there's a process they want to delete
    OptionSelected = False
    while OptionSelected == False:
        DeleteProcesses = input("Delete Processes? (y/n): ").lower().replace(' ','')
        if DeleteProcesses == 'y' or DeleteProcesses == 'n':
            OptionSelected = True
        else:
            print('\nInvalid Option\n')

    # If the user selected 'n', the connection is closed and the program terminates. Otherwise, the user
    # is prompted for the ID they wish to terminate
    if DeleteProcesses == 'y':
        # The program will ask until a valid ID is provided
        ProcessSelected = False
        while ProcessSelected == False:
            ProcessToDelete = input('\nEnter Process ID: ')
            # Just in case the user needs to shut things down, they can use the key 'q' to exit the program
            if ProcessToDelete == 'q':
                sys.exit(0)
            # IDs are integers and are stored as such in the process output, so the program makes sure the
            # input can be converted into an INT
            try:
                ProcessToDelete = int(ProcessToDelete)
            except ValueError:
                print('Process ID Not an Integer.')
            # The program also makes sure the process exists
            if ProcessToDelete in ProcessIDs:
                ProcessSelected = True
        # Once a valid ID has been selected, the program kills it
        mycursor.execute("KILL QUERY %s"%ProcessToDelete)
cnx.close()
                
         
