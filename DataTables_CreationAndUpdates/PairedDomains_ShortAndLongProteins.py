import os, sys, mysql.connector, datetime, random


'''
Authors : Jennifer James and Sara Willis
Date    : June 7, 2019

The purpose of this script is to collect tuples of Pfams and their metrics based on the lengths of their parent proteins. Two versions of each Pfam in the database are collected:

   - An instance of a pfam found in the shortest protein that contains it
   - An instance of a pfam found in the longest protein that contains it

If there are multiple instances of a pfam found in a protein/proteins of the same length which happen to be the longest/shortest, one is chosen at random. This can either happen if two proteins share the same length and both contain the pfam, if multiple instances of a pfam show up in a single protein, or some combination of the two.
'''



##########################################################################################
#                                   User-Specific Data                                   #
##########################################################################################

# User-Specific MySQL Connection Information
Database = ''
User     = ''
Host     = ''
Password = ''


# Unique Proteins Table
# -- More than one table name can be specified as a comma-delimited set of strings
# It's important, however, that the tables listed share the same structure and should
# be presented in the same order as their matching counterparts under Domain Metrics
# Table and Protein Metrics Table
# The full data tables probelmatically contain multiple transcripts per protein. Work
# was done by Catherine Weibel to select only one transcript per protein and this has
# been stored in the unique protein tables. 
UniqueProteinTables               = 'UniqueProteinGenePairs_Table_Ensembl','UniqueProteinGenePairs_Table_NCBI'
UniqueProteinUID                  = 'UID'

# Domain Metrics Table
DomainTables                      = 'Genomes_Multicellular_DomainMetrics','NCBIGenomes_DomainMetrics_Complete'
DomainProteinCorrespondance       = 'ProteinTableUID'
DomainTableUID                    = 'UID'
DomainISD                         = 'MeanISD_IUPred2_WithCys'
DomainAgg                         = 'DensityOfAggProneRegions'
DomainClustering                  = 'NormalizedIndexOfDispersion_Trunc_FILVM'
DomainLength                      = 'DomainLength'
PfamUID                           = 'PfamUID'

# Protein Metrics Table
ProteinMetricsTables              = 'Genomes_Multicellular_ProteinMetrics','NCBIGenomes_ProteinMetrics_Complete'
ProteinLength                     = 'ProteinLength'
ProteinTable_MetricCorrespondance = 'ProteinTableUID'

# Upload Table -- Where 
LengthContrastTable               = "PfamLengthContrastPairs"
Contrast_PfamUID                  = 'PfamUID'
Contrast_ShortProteinLength       = 'ShortProteinLength'
Contrast_ShortDomainUID           = 'ShortUID'
Contrast_ShortProteinTableUID     = 'ShortProteinTableUID'
Contrast_LongProteinLength        = 'LongProteinLength'
Contrast_LongDomainUID            = 'LongUID'
Contrast_LongProteinTableUID      = 'LongProteinTableUID'
Contrast_LongISD                  = 'LongISD'
Contrast_ShortISD                 = 'ShortISD'
Contrast_LongClustering           = 'LongClustering'
Contrast_ShortClustering          = 'ShortClustering'
Contrast_LongDensityAPRs          = 'LongDensityAPRs'
Contrast_ShortDensityAPRs         = 'ShortDensityAPRs'





##########################################################################################
#                                Program Executes Below                                  #
##########################################################################################

start_time = datetime.datetime.now()
print('Beginning Analysis\nCurrent Time: %s\n'%start_time)
sys.stdout.flush()

# A connection is established to the SQL database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# Pfam entries are stored in a nested dictionary. There will be entries for both the longest
# and shortest containing genes
PfamLongAndShortDictionary = {}

# This script is set up so that multiple tables can be pulled and their data are concatenated
# together. It's also possible to only pull from one table. The script makes sure the table
# names are in the correct format before the script is run. 
if type(UniqueProteinTables) == tuple:
    DataTables = zip(UniqueProteinTables,DomainTables,ProteinMetricsTables)
else:
    DataTables = [UniqueProteinTables,DomainTables,ProteinMetricsTables]


for DataTableSet in DataTables:
    UniqueProteinTable,DomainTable,ProteinMetricsTable = DataTableSet
    
    # A designator, e.g. 'Ensembl' or 'NCBI' is pulled from the table so that UIDs associated
    # with domains and their proteins can be distinguished from one another in the final output.
    # This allows proteins to be traced to their originating table
    Source = UniqueProteinTable.split('_')[-1]
    print('Extracting Entries from %s Data Tables\n'%Source)
    sys.stdout.flush()

    # The relevant tables are extracted using join statements. The unique table ensures we don't
    # get entries from multiple transcripts of the same gene (i.e. utilizes Catherine's work)
    ExtractionStatement = "SELECT dom.%s"%(',dom.'.join([DomainTableUID,DomainISD,DomainAgg,DomainClustering,DomainLength,PfamUID]))+",pmet.%s"%',pmet.'.join([ProteinLength,ProteinTable_MetricCorrespondance]) + " FROM "+ProteinMetricsTable+" AS pmet INNER JOIN " + UniqueProteinTable+ " AS u ON u.%s=pmet.%s"%(UniqueProteinUID,ProteinTable_MetricCorrespondance)+ " INNER JOIN "+DomainTable+" AS dom ON pmet.%s=dom.%s"%(ProteinTable_MetricCorrespondance,DomainProteinCorrespondance)
    
    intermediate_time = datetime.datetime.now()
    
    mycursor.execute(ExtractionStatement)
    results = mycursor.fetchall()
    
    print('%s Entries Successfully Extracted\nTime Taken: %s\nSorting Extracted Data\n'%(len(results),datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
    
    intermediate_time = datetime.datetime.now()

    for result in results:
        pfam_tableuid,pfam_isd,pfam_agg,pfam_clustering,domain_length,pfam_uid,protein_length,protein_uid = result
        
        # The UIDs from each table are transformed so their source table is obvious
        pfam_tableuid = str(pfam_tableuid)+ '_%s'%Source
        protein_uid = str(protein_uid)+ '_%s'%Source

        # And a similar tuple is designed to look like the extracted data with updated
        # identifying UIDs
        identifyingResult = (pfam_tableuid,pfam_isd,pfam_agg,pfam_clustering,domain_length,pfam_uid,protein_length,protein_uid)

        # If the pfam is not already in our output dictionary we add it as both the shortest
        # and the longest. Tuples will be saved as a list so that in the event of multiple
        # domains originating from a protein/proteins classified either as the shortest or
        # longest, one may be chosen at random
        if pfam_uid not in PfamLongAndShortDictionary:

            PfamLongAndShortDictionary[pfam_uid] = {'short': {'length' : protein_length,
                                                              'entry' : [identifyingResult]},
                                                    'long' : {'length' : protein_length,
                                                              'entry' : [identifyingResult]}}

        else:
            # If a protein is found to have a length that is greater than (less than) the
            # length of the current longest (shortest) protein containing the current pfam,
            # it replaces the previous entry
            if protein_length > PfamLongAndShortDictionary[pfam_uid]['long']['length']:
                PfamLongAndShortDictionary[pfam_uid]['long']['length'] = protein_length
                PfamLongAndShortDictionary[pfam_uid]['long']['entry']  = [identifyingResult]
            if protein_length < PfamLongAndShortDictionary[pfam_uid]['short']['length']:
                PfamLongAndShortDictionary[pfam_uid]['short']['length'] = protein_length
                PfamLongAndShortDictionary[pfam_uid]['short']['entry'] = [identifyingResult]
                
            # If a protein shares the same length as the current longest (shortest), the
            # result is added to the list 
            if protein_length == PfamLongAndShortDictionary[pfam_uid]['long']['length']:
                PfamLongAndShortDictionary[pfam_uid]['long']['entry'].append(identifyingResult)
            if protein_length == PfamLongAndShortDictionary[pfam_uid]['short']['length']:
                PfamLongAndShortDictionary[pfam_uid]['short']['entry'].append(identifyingResult)
    print('Entries Sorted\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
intermediate_time = datetime.datetime.now()

print('All Entries Successfully Sorted\nUploading Data\n')
sys.stdout.flush()

# Once all the entries have been accounted for, they are uploaded to the relevant MySQL table.
# Each result corresponds with a particular Pfam and the shortest/longest results are parsed
for entry in PfamLongAndShortDictionary:

    # Short Results <-- one entry is selected at random from the list of results
    short_pfam_tableuid,short_pfam_isd,short_pfam_agg,short_pfam_clustering,short_domain_length,pfam_uid,short_protein_length,short_protein_uid = random.choice(PfamLongAndShortDictionary[entry]['short']['entry'])

    # Long Results <-- the same applies here
    long_pfam_tableuid,long_pfam_isd,long_pfam_agg,long_pfam_clustering,long_domain_length,pfam_uid,long_protein_length,long_protein_uid= random.choice(PfamLongAndShortDictionary[entry]['long']['entry'])
    
    # All columns in the results table will be stored so a join statement can be used to
    # correctly format them for a SQL insertion statement
    InsertionColumns = [Contrast_PfamUID,Contrast_ShortProteinLength,Contrast_ShortDomainUID,Contrast_ShortProteinTableUID,Contrast_LongProteinLength,Contrast_LongDomainUID,Contrast_LongProteinTableUID,Contrast_LongISD,Contrast_ShortISD,Contrast_LongClustering,Contrast_ShortClustering,Contrast_LongDensityAPRs,Contrast_ShortDensityAPRs]
    InsertionStatement = "INSERT INTO "+LengthContrastTable+" ("+','.join(InsertionColumns)+") VALUES ("+','.join(['%s' for i in InsertionColumns]) +")"
    InsertionData = (pfam_uid,short_protein_length,short_pfam_tableuid,short_protein_uid,long_protein_length,long_pfam_tableuid,long_protein_uid,long_pfam_isd,short_pfam_isd,long_pfam_clustering,short_pfam_clustering,long_pfam_agg,short_pfam_agg)
    mycursor.execute(InsertionStatement,InsertionData)
    cnx.commit()

print('Entries Successfully Uploaded\nTime Taken: %s\n\nScript Complete!\nTotal Time Taken: %s'%(datetime.datetime.now()-intermediate_time,datetime.datetime.now()-start_time))
sys.stdout.flush()

cnx.close()
