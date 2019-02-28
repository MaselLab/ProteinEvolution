import os, sys, json, csv, mysql.connector, datetime, ftputil
from Bio.Seq import Seq
from Bio import SeqIO


'''
Author : Sara Willis
Date   : February 11, 2019
--------------------------


This script is used to access the Ensembl databases so that it can extract IDs associated with chloroplast and mitochondrial genes.

This script was written so that mitochondrial and chloroplast genes could be identified in our dataset containing the coding sequences from over 400 genomes. This is in part because mitochondrial and chloroplast genes tend to use alternative coding alphabets. 

This script connects to the Ensembl FTP sites, checks various species names located in those databases and cross-references them with our databases. In the case that the species in the FTP site is one in our database, we extract that species, otherwise it's ignored. 

Once a species' Fasta file containing its coding sequences is downloaded, the description of each sequence is parsed. What we're looking for is either the identifier 'Pt' or 'MT' under 'chromosome:' that identifies the sequence as either chloroplast or mitochrondrial, respectively. 

If a chloroplast or mitochrondrial sequence is found, its UID is recorded in an output file for later use.
'''

################################################################################
#                          User-Specific Information                           #
################################################################################


# User-supplied MySQL connection information
Database = ''
User = ''
Host = ''
Password = ''

# MySQL data table where species data are stored
SpeciesDatabase = 'SpeciesList'

# The column name where the species' names are stored in Newick format
SpeciesNameColumn = 'NewickSpeciesName'

# The column that describes where the species data came from (i.e. V93)
SourceColumn = 'SourceDatabase'

# When Verbose is set to True, the program will keep the user updated on the progress of the script
Verbose = True


################################################################################
#                           Program Executes Below                             #
################################################################################


SpeciesList = {}

start_time = datetime.datetime.now()

# The filenames that will be saved for the user
ChloroplastOutputFile = open('EnsemblChloroplastGenes.txt','w')
MitochondrialOutputFile = open('EnsemblMitochondrialGenes.txt','w')

# A connection to the user's MySQL database is established
if Verbose == True:
    print('\nConnecting to MySQL database\n\n')
    sys.stdout.flush()
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
# And a cursor is defined so we can interact with the data tables
mycursor = cnx.cursor(buffered = True)

# We want the names of the species that are located in our database
pullSpeciesStatement = "SELECT "+SpeciesNameColumn+","+SourceColumn+" FROM " + SpeciesDatabase
mycursor.execute(pullSpeciesStatement)

# We extract the species names (in Newick format, which is how genomes are stored in the FTP repositories)
# in our database, as well as where they were extracted from. This helps us pinpoint where to search for specific genomes
for element in mycursor:
    SpeciesName = element[0].lower()
    Database = element[1]
    if Database not in SpeciesList:
        SpeciesList[Database] = [SpeciesName]
    else:
        SpeciesList[Database].append(SpeciesName)
if Verbose == True:
    print('Species data extracted\nTime Taken: %s\n\nConnecting to Ensembl V93 FTP repository'%(datetime.datetime.now()-start_time))
    sys.stdout.flush()
intermediate_time = datetime.datetime.now()

# We first start with Ensembl Version 93 by connecting to Ensembl's FTP site
FTPHost = ftputil.FTPHost('ftp.ensembl.org','anonymous','password')
FTPHost.chdir('/pub/release-93/fasta/')

# The list of directories gives us the names of the species so comparisons can be made
# with our species, and also so we know the names of the folders we need to navigate into
dir_list = FTPHost.listdir(FTPHost.curdir)
for directory in dir_list:
    
    # If the name of the species (which is the directory listing) is in our database, we extract
    # the coding sequences from that species (directory)
    if directory in SpeciesList['V93']:
        FTPHost.chdir('//pub/release-93/fasta/%s/cdna/'%directory)
        
        # We don't know specifically what the names of the files are in the directory, so we go
        # through everything that's located in our current directory and look for the identifier
        # 'cdna.all.fa.gz' which is common to all coding sequence fasta filenames
        files = FTPHost.listdir(FTPHost.curdir)
        targetFile = [i for i in files if 'cdna.all.fa.gz' in i][0]
        
        # We only try downloading the coding sequence file if it exists
        if targetFile != '':
            
            # If a previous run has downloaded and unzipped the file, we don't need to download it again
            # so we skip the download/unzipping process
            if os.path.exists(targetFile.replace('.gz','')) == False:
                if Verbose == True:
                    print('Downloading %s from EnsemblV93'%targetFile)
                    sys.stdout.flush()
                FTPHost.download(targetFile,targetFile)
                os.system('gzip -d %s' %targetFile)
                
            # The unzipped file no longer has the .gz extension
            fastaFilename = targetFile.replace('.gz','')

            # Each entry in the fasta file is then searched for either chloroplast or mitochondrial genes
            # using the identifiers MT (mitochondrial) or Pt (chloroplast). If one is found, it's written
            # into the output file
            for record in SeqIO.parse(fastaFilename,'fasta'):
                description = record.description
                if 'chromosome:' in description:
                    chromosomeStart = str(description).index('chromosome:')
                    chromosome = description[chromosomeStart:].split(' ')[0].replace('chromosome:','')
                    if 'Pt' in chromosome:
                        ChloroplastOutputFile.write('%s\n'%record.id.split('.')[0])
                    elif 'MT' in chromosome:
                        MitochondrialOutputFile.write('%s\n'%record.id.split('.')[0])
                        
            # We want to mitigate clutter, so we remove fasta files once we're done with them
            os.remove(fastaFilename)
if Verbose == True:
    print('\n\nMitochondrial and Chloroplast genes successfully extracted from Ensembl V93 FTP repository\n\nTime taken: %s'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()


# We repeat the same process above with the repositories for V40 Plants, Protists, Metazoa, and Fungi
targetDirectories = [('fungi','EnsemblFungiV40'),('protists','EnsemblProtistV40'),('metazoa','EnsemblMetazoaV40'),('plants','EnsemblPlantV40')]
if Verbose == True:
    print('Connecting to V40 Ensembl repositories\n')
    sys.stdout.flush()
FTPHost = ftputil.FTPHost('ftp.ensemblgenomes.org','anonymous','password')
for targetDirectory in targetDirectories:
    intermediate_time = datetime.datetime.now()
    FTPHost.chdir('//pub/%s/release-40/fasta/'%targetDirectory[0])
    dir_list = FTPHost.listdir(FTPHost.curdir)
    for directory in dir_list:
        if directory in SpeciesList[targetDirectory[1]]:
            FTPHost.chdir('//pub/%s/release-40/fasta/%s/cdna/' %(targetDirectory[0],directory))
            files = FTPHost.listdir(FTPHost.curdir)
            targetFile = [i for i in files if 'cdna.all.fa.gz' in i][0]
            if Verbose == True:
                print('Downloading %s from repository %s'%(targetFile,targetDirectory[1]))
                sys.stdout.flush()
            if targetFile != '':
                if os.path.exists(targetFile.replace('.gz','')) == False:
                    FTPHost.download(targetFile,targetFile)
                    os.system('gzip -d %s'%targetFile)
                fastaFilename = targetFile.replace('.gz','')
                for record in SeqIO.parse(fastaFilename,'fasta'):
                    description = record.description
                    if 'chromosome:' in description:
                        chromosome = description[description.index('chromosome:'):].split(' ')[0].replace('chromosome:','')
                        if 'Pt' in chromosome:
                            ChloroplastOutputFile.write('%s\n'%record.id.split('.')[0])
                        elif 'MT' in chromosome:
                            MitochondrialOutputFile.write('%s\n'%record.id.split('.')[0])
                            
                            
                os.remove(fastaFilename)
    if Verbose == True:
        print('Genes successfully parsed from %s\nTime taken: %s\n\n'%(targetDirectory[0],datetime.datetime.now()-intermediate_time))
        sys.stdout.flush()
        
          
ChloroplastOutputFile.close()
MitochondrialOutputFile.close()

if Verbose == True:
    print('Mitochondrial and chloroplast files successfully created\nTotal time taken: %s\n\n'%(datetime.datetime.now()-start_time))
    sys.stdout.flush()

