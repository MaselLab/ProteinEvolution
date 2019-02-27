from Bio import SeqIO
from Bio.Seq import Seq
import os, sys, json, csv, re, string, ftputil, mysql.connector, datetime


OutputAccessionFile = open('NCBI_mitochondrial_and_chloroplast_geneIDs.txt','w')

'''

Author : Sara Willis
Date   : Wednesday February 27, 2019


The purpose of this script is to search for mitochondrial or chloroplast genes in full genome fasta files on the NCBI FTP site. Only species in our dataset will be searched and all accession numbers associated with the chloroplast/mitochrondrial genes will be written to a file
'''

########################################################################
######################### User-Specific Data    ########################
########################################################################

Database = ''
User = ''
Host = ''
Password = ''

# Temporary Directory
TempFilesDirectoryPath = './TempFiles/'
TempFilesDirectoryFilename = TempFilesDirectoryPath.split('/')[1]

# Temporary Fasta Files
ReformattedDownloadedFile = TempFilesDirectoryPath + 'Reformatted.fa'
TempFastaFile = TempFilesDirectoryPath + 'Temp.fa'
TempIUPredFastaFile = 'Temp.fa'

Verbose = True


########################################################################
############################ Submodules ################################
########################################################################

def AttemptFileDownload(directory, subdirectory, genomeType, TempFilesDirectoryPath):
    host = ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')
    # Change to the indicated genomeType
    host.chdir('/genomes/refseq/%s/%s/%s'%(directory, subdirectory, genomeType))
    if Verbose == True:
        print('\nDirectory Found')
        sys.stdout.flush()
    # Determinet the name of the one included subdirectory
    contentOfSpeciesDirectory = host.listdir(host.curdir)[0]
    # navigate into that subdirectory
    host.chdir('/genomes/refseq/%s/%s/%s/%s'%(directory, subdirectory,genomeType,contentOfSpeciesDirectory))
    # determine if the directory has coding sequences associated with it
    filesList = host.listdir(host.curdir)
    FeaturesFileExists = False
    for File in filesList:
        # Checks each file/directory in the current location for a file
        # with the extension that indicates the coding sequence file
        if '_feature_table.txt.gz' in File:
            FeaturesFile = File
            # If the file found, the loop is exited
            FeaturesFileExists = True
            break
        # If the coding sequence file exists, it is downloaded
    if FeaturesFileExists == True:
        if Verbose == True:
            print('\nDownloading File...\n')
            sys.stdout.flush()
            
        # The file is saved in the temporary files directory
        TargetFilename = TempFilesDirectoryPath + FeaturesFile
        # Downloads the file
        host.download(FeaturesFile, TargetFilename)
        # Checks to see that the file was downloaded
        if os.path.exists(TargetFilename) == True:
            if Verbose == True:
                print('\nFile Successfully Downloaded.\n')
                print('\nUnzipping File...\n')
                sys.stdout.flush()
                
            # if the file was downloaded correctly, it is unzipped
            os.system('gzip -d %s' %TargetFilename)
            # Checks to see that the file was unzipped successfully
            UnzippedTargetFilename = TargetFilename.replace('.gz','')
            # Checks to make sure the file was unzipped successfully
            if os.path.exists(UnzippedTargetFilename) == True and Verbose == True:
                print('\nFile Successfully Unzipped.\n')
                sys.stdout.flush()
        else:
            print('There was a problem downloading the file')
        return UnzippedTargetFilename


 


########################################################################
##################### PROGRAM EXECUTES BELOW ###########################
########################################################################
start_time = datetime.datetime.now()
print('Beginning NCBI Extraction\n%s'%start_time)

SpeciesProcessed = 0

# Filename for relevant species
NCBISpeciesFilename =  'Species_NotInEnsembl_CompleteInGold_InTimetree_NotInSpeciesList.csv'

# Creates a list of the species that will attempt to be extracted
speciesList = {}

if os.path.isdir(TempFilesDirectoryFilename) == False:
    os.system('mkdir %s'%TempFilesDirectoryFilename)


cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
# Defines a cursor so we can interact with the datatable
mycursor = cnx.cursor(buffered=True)
mycursor.execute("SELECT SpeciesUID,NewickSpeciesName FROM SpeciesList WHERE SourceDatabase='NCBI'")
SpeciesResults = mycursor.fetchall()
for SpeciesResult in SpeciesResults:
    speciesList[SpeciesResult[1]] = SpeciesResult[0]
    
# A list of the directories included in the NCBI FTP site that will be searched
targetDirectories = ['vertebrate_other', 'vertebrate_mammalian', 'plant','invertebrate','fungi','protozoa']

# Connect to the NCBI FTP site
host = ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')

# Change to the genomes directory
host.chdir('/genomes/refseq/')

# List the directories inside 'genomes'
dir_list = host.listdir(host.curdir)

# Check each directory to see if it's located in the desired directories list
for directory in dir_list:
    # directory corresponds to vertebrate_other, vertebrate_mammalian, etc...
    if directory in targetDirectories:
        # If the directory is one we want to check, we navigate to it
        host = ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')

        # Change to target directory 
        #host.chdir('//genomes/refseq/%s'%directory)
        host.chdir('/genomes/refseq/%s'%directory)

        # We now check each subdirectory and compare it to our speciesList
        speciesDirectories = host.listdir(host.curdir)
        # subdirectories will have species names. Each subdirectory will be compared to our species list
        for speciesDirectory in speciesDirectories:            
            CodingSequencesDictionary = {}
            ProteinSequencesDictionary = {}
            # if the subdirectory is for a species we're looking forward, we
            # navigate into that subdirectory
            if speciesDirectory in speciesList:
                SpeciesUID = speciesList[speciesDirectory]

                downloadSuccessfullyExecuted = True
                # We try to navigate to the representative genome
                try:
                    if Verbose == True:
                        print('Attempting to Download Representative Genome')
                        sys.stdout.flush()
                    UnzippedTargetFilename = AttemptFileDownload(directory, speciesDirectory, 'representative', TempFilesDirectoryPath)
                    if Verbose == True:
                        print(UnzippedTargetFilename)
                        sys.stdout.flush()
                except:
                    if Verbose == True:
                        print('Representative Genome Not Found\nAttempting to Download Latest Genome Version')
                        sys.stdout.flush()
                    try:
                        UnzippedTargetFilename = AttemptFileDownload(directory, speciesDirectory, 'latest_assembly_versions', TempFilesDirectoryPath)
                    except:
                        if Verbose == True:
                            print('No Coding Sequence Files Found\nSpecies Skipped')
                            sys.stdout.flush()
                        downloadSuccessfullyExecuted = False
                if downloadSuccessfullyExecuted == False:
                    pass
                else:
                    SpeciesProcessed += 1
                    # If the file was successfully acquired, it's parsed
                    Header = True
                    with open(UnzippedTargetFilename,'r') as f:
                        reader = csv.reader(f, delimiter = '\t')
                        for row in reader:

                            if Header == True:
                                Header = False
                            else:
                                Feature = row[0]
                                Class = row[1]
                                Assembly = row[2]
                                AssemblyUnit = row[3]
                                SeqType = row[4]
                                Chromosome = row[5]
                                GenomicAccession = row[6]
                                Start = row[7]
                                End = row[8]
                                Strand = row[9]
                                ProductAccession = row[10]
                                NonRedundantRefSeq = row[11]
                                RelatedAccession = row[12]
                                Name = row[13]
                                Symbol = row[14]
                                GeneID = row[15]
                                LocusTag = row[16]
                                FeatureIntervalLength = row[17]
                                ProductLength = row[18]
                                Attributes = row[19]
                                if ('MT' in Chromosome or 'Pltd' in Chromosome) and ProductAccession != '':
                                    OutputAccessionFile.write('%s\t%s\t%s\n'%(SeqType,Chromosome,ProductAccession))

                                
                                
                    os.remove(UnzippedTargetFilename)

                    
OutputAccessionFile.close()
print('Number of species processed: %s' %SpeciesProcessed)
print('Total time taken: %s' %(datetime.datetime.now()-start_time))
