from Bio import SeqIO
from Bio.Seq import Seq
import os, sys, json, csv, re, string, mysql.connector, datetime



'''
Author: Sara Willis
Date  : February 18, 2019

Note: I haven't added deletion statements. You may want to write a separate script to do this after this one has run successfully. This is because we don't want to accidentally delete data before we're *sure* we've transferred it successfully to the new table.

'''

# These are the data tables where the protist entries are stored
DataTables = ['EnsemblGenomes_Coding_Complete','EnsemblGenomes_DomainMetrics_Complete','EnsemblGenomes_Protein_Complete','EnsemblGenomes_ProteinMetrics_Complete']


# This logs you into MySQL, just put in your fusion access Username and password
Database = ''
User = ''
Host = ''
Password = ''

# This sets up a connection with MySQL using the information you provided above
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)

# a cursor allows you to interact with the database for extractions/uploads/updates etc. This is where all
# the SQL syntax is used
mycursor = cnx.cursor(buffered = True)

# We start by extracting the speciesUID, NewickSpeciesName, and Source Database from our list of species
mycursor.execute("SELECT SpeciesUID,NewickSpeciesName,SourceDatabase FROM SpeciesList_copy")

# fetchall() pulls the results out of our cursor as a list of tuples
SpeciesResults = mycursor.fetchall()

# Now we search through the species we've extracted
for speciesEntry in SpeciesResults:

    # We assign variable names for ease of use
    SpeciesUID,SpeciesName,Source = speciesEntry

    # We're only interested in protist genes, and since they're all in the Ensembl database, we identify them
    # easily since they all came from the ProtistV40 Ensembl repository. If we identify a species as a protist, we
    # set about moving those entries from our original databases to new protist-specific ones
    if Source == 'EnsemblProtistV40':

        # Genomic data is stored in four separate data tables, so we're going to need to move entries from all of them
        for datatable in DataTables:

            # It's a pain in the ass to write in all the column names, so we use this statement to extract the column
            # names from the source datatable. Since the tables we're moving our entries over to are structural copies
            # of the source tables, we don't need to worry about altering the column names for the sake of uploading
            mycursor.execute("SHOW COLUMNS FROM "+datatable)
            Columns = mycursor.fetchall()
            # Column names come with a lot of garbage, so we strip off all the unneccesary gunk
            Columns = [i[0] for i in Columns]

            # This just pulls which data table we're looking at (i.e. coding, protein, proteinmetrics, or domainmetrics
            # It uses this to upload to the appropriate protist data table
            datatableType = datatable.split('_')[1]
            InsertTable = "Genomes_Protists_%s"%datatableType

            # We define an instertion statement to get our entries into the protist data tables
            InsertionStatement = "INSERT INTO " +InsertTable+" ("+','.join(Columns)+") VALUES (" +','.join(['%s' for i in Columns])+")"

            # finally, we extract species-specific data from our source data table using the speciesUID
            ExtractionStatement = "SELECT * FROM "+datatable+" WHERE SpeciesUID=%s"%SpeciesUID
            mycursor.execute(ExtractionStatement)
            protistResults = mycursor.fetchall()
            
            # We then go through each row and upload it to the protist database
            for protistResults in protistResults:
                mycursor.execute(InsertionStatement,protistResults)
                cnx.commit()


cnx.close()
