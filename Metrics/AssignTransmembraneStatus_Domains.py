import os,sys, mysql.connector, datetime, re

'''
Authors : Sara Willis and Jennifer James
Date    : May 31, 2019

This script is designed to determine the transmembrane status of Pfam domains based on tmhmm predictions. 

Tmhmm predictions output strings of indices associated with helix boundaries. A domain is counted as transmembrane if either:

   1) At least 50% of the domain is contained with a helix

   or

   2) At least 50% of a helix is contained within a domain

These two conditions are necessary as only accepting condition 1) would reject long domains completely containing a helix. Similarly, only accepting codition 2) would reject short pfams completely contained within a long helix.

As the cutoff of 50% is arbitrary, this script does not assign binary values to the output but instead outputs a floating value which it uploads to the domain metric table allowing some fiddling in downstream analyses if desired.

The output value is the maximum fraction either of the ... and is calculated in the following way

   a) The boundaries of each helix are searched against the domain boundaries
   b) If overlap is found, the fraction of the overlap relative to both the domain and the helix is found and the max value is retained
   c) Once all overlaps are accounted for a particular domain, the maximum value of all retained fractions is selected and uploaded to the domain metrics table under TmhmmTopology

# This makes it so that downstream analysis is selecting a cutoff based on the best helix overlaps for each domain

'''

####################################################################################
#                                  User Options                                    #
####################################################################################

# User-specific MySQL connection information
Database = ''
User     = ''
Host     = ''
Password = ''

# The relevant tables where all data are stored. Multiple tables can be listed in comma-delimited format
# The linked tables should all appear in the same position for each of the tables otherwise the data
# will not be extracted correctly
ProteinTables       = 'NCBIGenomes_Protein_Complete'#,'Genomes_Multicellular_Protein'
ProteinMetricTables = 'NCBIGenomes_ProteinMetrics_Complete'#,'Genomes_Multicellular_ProteinMetrics'
DomainMetricTables  = 'NCBIGenomes_DomainMetrics_Complete'#,'Genomes_Multicellular_DomainMetrics'

# Column names where data are stored -- the column names should match in all tables listed above
PfamStartColumn               = 'PfamStart'           # Protein Table
PfamStopColumn                = 'PfamStop'            # Protein Table
PfamUIDsColumn                = 'PfamMetricsTableUID' # Protein Table
TopologyColumn                = 'Topology'            # Protein Metrics Table
ProteinTableUIDColumn         = 'UID'                 # Protein Table
ProteinMetricProteinUIDColumn = 'ProteinTableUID'     # Protein Metrics Table -- links metrics and protein table
DomainTransmembraneColumn     = 'TmhmmTopology'       # Domain Metrics Table  -- Where output will be stored
DomainUIDColumn               = 'UID'                 # Domain Metrics Table



####################################################################################
#                             Program Executes Below                               #
####################################################################################

start_time = datetime.datetime.now()
print('\nProgram Executing\nCalculating TmhmmTopology of Pfam Domains\nCurrent Time: %s\n'%start_time)
sys.stdout.flush()

# to allow for the possibility that we're only updating one table instead of multiple,
# we check to make sure our data are in the correct format
if type(ProteinTables) == tuple:
    # The tables are combined so they can be easily used together
    ZippedTables = zip(ProteinTables,ProteinMetricTables,DomainMetricTables)
# If not, we force them to be so the program can correctly execute
else:
    ZippedTables = [(ProteinTables,ProteinMetricTables,DomainMetricTables)]


# A connection is established to the MySQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Each database is then analyzed independently
for Tables in ZippedTables:
    intermediate_time = datetime.datetime.now()
    ProteinTable,ProteinMetricsTable,DomainTable = Tables
    print('Extracting Data from %s'%','.join(Tables))
    # The relevant information to calculate the transmembrane status of each domain in
    # the database is pulled from the MySQL tables
    ExtractionStatement = "SELECT protein."+PfamStartColumn+",protein."+PfamStopColumn+",protein."+PfamUIDsColumn+",pmetric."+TopologyColumn+" FROM "+ProteinTable+" AS protein INNER JOIN "+ProteinMetricsTable+" AS pmetric ON protein."+ProteinTableUIDColumn+"=pmetric."+ProteinMetricProteinUIDColumn#+" WHERE protein.UID < 1000"

    intermediate_time=datetime.datetime.now()

    mycursor.execute(ExtractionStatement)
    print('Data Successfully Extracted\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
    results = mycursor.fetchall()

    # Each entry is specifically for a full  protein so nested loops are needed to analyze
    # the associated domains
    print('Calculating Transmembrane Status of Domains')
    sys.stdout.flush()
    
    for result in results:  
        PfamStarts,PfamStops,PfamUIDs,Topology = result
        # The boundaries of transmembrane helices are output from Tmhmm as strings separated by '-'
        # and differentiated by either 'i' or 'o'. The indices of the boundaries are pulled, grouped 
        # and stored as integers within sublists in a master list
        transmembrane_helices = [list(map(int,j)) for j in [i.split('-') for i in re.findall('[0-9]+-[0-9]+', Topology)] ]#for k in j]
        # Pfam indices are stored in the protein tables as comma-delimited strings. They are stored in
        # a similar way to the transmembrane helices
        PfamIndices = [list(map(int,j)) for j in zip(PfamStarts.split(','),PfamStops.split(','),PfamUIDs.split(','))]

        # Each domain is then searched against the transmembrane boundaries
        for index in PfamIndices:

            # The max fraction of each pfam/helix overlap is saved to a list for later analysis
            HelixOverlaps = []

            # To get the length of the overlap, we use set intersections
            PfamSpan = set(range(index[0],index[1]))

            # If there are no transmembrane helices, then the overlap is automatically zero. Otherwise,
            # we have to check each helix against the domain boundaries
            if len(transmembrane_helices) != 0:
                for helix in transmembrane_helices:
                    HelixSpan = set(range(helix[0],helix[1]))
                    overlap = PfamSpan.intersection(HelixSpan)
                    if len(overlap) != 0:
                        FractionPfam = len(overlap)/len(PfamSpan)
                        FractionHelix = len(overlap)/len(HelixSpan)
                        MaxPercentage = max(FractionPfam,FractionHelix)
                        HelixOverlaps.append(MaxPercentage)
                    else:
                        HelixOverlaps.append(0)
            else:
                HelixOverlaps.append(0)

            # Once all the fractional values have been determined, we pull the max and store it in the domain
            # metrics table
            MaxCoverage = max(HelixOverlaps)
            UpdateStatement = "UPDATE " + DomainTable + " SET "+DomainTransmembraneColumn+"=%s WHERE "+DomainUIDColumn+"=%s"
            UpdateData = (MaxCoverage,index[2])
            mycursor.execute(UpdateStatement,UpdateData)
            cnx.commit()
    print('Transmembrane Status Successfully Completed\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()

print('Program Complete\nTime Taken: %s\nCurrent Time: %s\n'%(datetime.datetime.now()-start_time,datetime.datetime.now()))
sys.stdout.flush()
cnx.close()
