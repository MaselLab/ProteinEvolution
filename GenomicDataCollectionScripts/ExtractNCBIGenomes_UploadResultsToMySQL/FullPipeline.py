from DataExtractionSubmodules import RemoveAllFiles, CreateFastaFile, CreateNCBISpeciesList, RunInterProScan, CalculateIUPred, LoadConnectionInformation, CodingSequenceQualityControl, AttemptFileDownload, FormatAndExtractSpeciesUID
from Bio import SeqIO
from Bio.Seq import Seq
import os, sys, json, csv, re, string, ftputil, mysql.connector, datetime

#There is a bug somewhere in this code that uploads multiple copies of a single protein. This needs to be fixed

'''
This script does the following:

   1) Extracts coding sequences from the NCBI databases
     - Coding sequences are for species that pass through the following filters:
       a) Species not in Ensembl Databases
       b) Species listed as fully-sequenced in GOLD database
       c) Species found in TimeTree.org
       d) Metazoa and Viridiplantae only
     - All species come from a supplied CSV file

   2) Converts coding sequences to protein sequences

   3) Runs them through InterProScan to get the Pfams associated with each gene
     - Genes with no Pfam hits are discarded from analyses

   4) Runs sequences with Pfam hits through IUPred to get disorder predictions
     - Removes cysteines from all protein sequences before analysis
     - Calculates average ISD for full protein sequence
     - Calculates average ISD for all Pfam domains associated with the sequence
     - Stores raw IUPred output

   5) Uploads data to a MySQL database
'''

########################################################################
######################### Load MySQL Login Data ########################
########################################################################

Database, User, Host, Password = LoadConnectionInformation()

########################################################################


########################################################################
###################### Temporary Paths and Filenames ###################
########################################################################

# Temporary Directory
TempFilesDirectoryPath = './TempFiles/'
TempFilesDirectoryFilename = TempFilesDirectoryPath.split('/')[1]

# Removes Temporary Directory and contents and creates empty directory
# before every run
RemoveAllFiles(TempFilesDirectoryPath)

# Temporary InterProScan filenames
InterProScanExecutablePath = '../bin/interproscan-5.30-69.0/interproscan.sh'
InterProScanOutputPath = TempFilesDirectoryPath + './InterProScanOutput.txt'

# Temporary Fasta Files
ReformattedDownloadedFile = TempFilesDirectoryPath + 'Reformatted.fa'
TempFastaFile = TempFilesDirectoryPath + 'Temp.fa'
TempIUPredFastaFile = 'Temp.fa'


########################################################################
##################### PROGRAM EXECUTES BELOW ###########################
########################################################################

SpeciesProcessed = 0

# Filename for relevant species
NCBISpeciesFilename =  'Species_NotInEnsembl_CompleteInGold_InTimetree_NotInSpeciesList.csv'

# Creates a list of the species that will attempt to be extracted
speciesList = CreateNCBISpeciesList(NCBISpeciesFilename)

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
                SpeciesUID = FormatAndExtractSpeciesUID(speciesDirectory)

                downloadSuccessfullyExecuted = True
                # We try to navigate to the representative genome
                try:
                    print('Attempting to Download Representative Genome')
                    UnzippedTargetFilename = AttemptFileDownload(directory, speciesDirectory, 'representative', TempFilesDirectoryPath)
                    print(UnzippedTargetFilename)
                except:
                    # If the representative genome is not available, then the latest assembly version is searched for
                    print('Representative Genome Not Found\nAttempting to Download Latest Genome Version')
                    try:
                        UnzippedTargetFilename = AttemptFileDownload(directory, speciesDirectory, 'latest_assembly_versions', TempFilesDirectoryPath)
                    except:
                        # If the latest assembly version also doesn't exist, then the species is skipped
                        print('No Coding Sequence Files Found\nSpecies Skipped')
                        downloadSuccessfullyExecuted = False
                if downloadSuccessfullyExecuted == False:
                    pass
                else:
                    SpeciesProcessed += 1
                    # If the file was successfully acquired, it's parsed
                    n= 0 
                    for record in SeqIO.parse(UnzippedTargetFilename, 'fasta'):
                        # The Gene ID and Protein ID are extracted from the gene description if it's
                        # possible to do so (some files don't contain this information). If either doesn't exist
                        # then they are just defined as empty strings
                        try:
                            GeneID = record.description[record.description.index('GeneID'):].split(']')[0].split(':')[1]
                        except:
                            GeneID = ''
                        try:
                            ProteinID = record.description[record.description.index('protein_id'):].split(']')[0].split('=')[1]
                        except:
                            ProteinID = ''
                        # The record accession is extracted
                        SequenceAccession = record.id
                        # The coding sequence is extracted and converted from a Sequence object to a string object
                        CodingSequence = str(record.seq)
                        # The protein sequence is found by converting the coding sequence. The stop codon(s) is(are) removed
                        ProteinSequence = str(record.seq.translate()).replace('*','')
                        # A unique identifier is made by combining the protein, gene and accession IDs
                        # The ID will be used to store relevant information about each sequence is different dictionaries
                        UID = GeneID + '||' + ProteinID + '||' + SequenceAccession
                        CodingSequencesDictionary[UID] = CodingSequence
                        ProteinSequencesDictionary[UID] = ProteinSequence

                        # @@@@@@@@@@ FOR TESTING SMALL NUMBER OF SEQUENCES @@@@@@@@@ #
                        # To ensure the pipeline is functioning before running it on the full genomes of all species,
                        # this allows the user to test the pipeline on a small number of sequences
                        if n >= 10:
                            break
                        # A simple fasta file is made to replace the original downloaded document
                        # This is done to replace the coding sequences with protein sequences, and replaces the
                        # original IDs with the new ones so that when the sequences are run through
                        # InterProScan, the results are easier to parse
                        CreateFastaFile(UID, ProteinSequence, ReformattedDownloadedFile,speciesDirectory)
                        n+= 1
                    # Once the new fasta file is made, the original coding sequence file is deleted
                    os.remove(UnzippedTargetFilename)

                    '''
                    The fasta file is then run through InterProScan to check each protein sequence against Pfam IDs. This exactly replicates how Pfam IDs were assigned to the genes in Ensembl. The results of the scan are saved as a dictionary. 

                    The dictionary of the InterProScan results is returned to the user as output from the function. Each entry is stored with the UID as the key, with a nested dectionary as the entry. Each nested dictionary has four entries: a list of Pfam UIDs associated with the protein encoded by that gene, the starting indices of that Pfam domain in the body of the protein sequence, the stopping indices, and the Evals associated with each Pfam hit. Each Pfam UID and its respective start/stop indices is stored with the same index in their respective lists so that they can be grouped and easily parsed.

                    Generic form of output dictionary:

                    {GeneUID: {Pfams: [Pfam_1, ..., Pfam_n], PfamStart: [Pfam_1_start, ..., Pfam_n_start], PfamStop = [Pfam_1_stop, ..., Pfam_n_stop], Eval = [Pfam_1_eval, ..., Pfam_n_eval]}}

                    Each GeneUID is defined above as the string GeneID, ProteinID, SequenceAccession joined by the delimiter '||'
                    This same key will be used in both the coding sequence and protein sequence dictionaries so all information can be joined together for the final output. 
                    '''

                    # InterProScan is run and the results are saved to a dictionary
                    InterProScanResultsDictionary = RunInterProScan(ReformattedDownloadedFile, InterProScanExecutablePath,InterProScanOutputPath)

                    # After InterProScan is run on the protein sequences, the Fasta file that was analyzed is removed
                    if os.path.exists(ReformattedDownloadedFile) == True:
                        os.remove(ReformattedDownloadedFile)

                    '''
                    The next major block of code uses the PfamUIDs, Coding Sequence Dictionary and Protein Sequence dictionary to calculate the ISD of the full sequences and the Pfam regions using IUPred (calculated using the 'long' option). 

                    Only sequences that have a PfamUID associated with them are used. Sequences without PfamUIDs assigned to them are discarded from the dataset. 

                    While all sequences with Pfam UIDs will be kept for phylogenetic purposes (i.e. the ages of each PfamUID will be dated using the entire set of sequences and their associated PfamUIDs, a stricter filter is implemented for calculating ISD. Namely, sequences that do not start with a start codon, have in-frame stop codons, or are not a multiple of three in length are discarded from the ISD calculations. This is done using a function that filters each sequence called CodingSequenceQualityControl(). The function takes in a coding sequence to analyze and  returns either True or False based on whether the sequence passed the required filters. If it returns True, the protein sequence associated with that coding sequence is written to its own fasta file and is fed into IUPred. If False is returned (i.e. it didn't pass the quality checks), it is skipped.
                    '''

                    # The time is printed for the user for when IUPred starts to run
                    currentTime = datetime.datetime.now()
                    print('\n%s: Starting IUPred Analysis\n' %currentTime)
                    # Only sequences that have an associated Pfam UID associated with them get passed on to be analyzed/uploaded
                    # to our database
                    for CodingID in InterProScanResultsDictionary:
                        # The coding sequence is pulled from the relevant dictionary
                        CodingSequence = CodingSequencesDictionary[CodingID]
                        # A quality control is implemented to check whether a coding sequence meets the following criteria:
                        #   1) Starts with a start codon
                        #   2) Has no in-frame stop codons
                        #   3) Is a multiple of three in length
                        SequencePassed = CodingSequenceQualityControl(CodingSequence)
                        # If the sequence meets the requirements, it is run through IUPred
                        if SequencePassed == True:
                            # First, the relevant information is pulled from various dictionaries
                            # The Protein sequence is located
                            ProteinSequence = ProteinSequencesDictionary[CodingID]
                            # The PfamUID is found
                            PfamIDs = InterProScanResultsDictionary[CodingID]['PfamUID']
                            # The starting and stopping locations for each Pfam UID are pulled
                            PfamStarts = InterProScanResultsDictionary[CodingID]['PfamStart']
                            PfamStops = InterProScanResultsDictionary[CodingID]['PfamStop']
                            # And the Evals for each Pfam UID are also pulled
                            InterProScanEval = InterProScanResultsDictionary[CodingID]['Eval']
                            # Next, IUPred is run on each protein sequence using the CalculateIUPred submodule
                            # The user is returned the full average ISD, as well as the average ISD for each
                            # Pfam region, as well as the full list of raw scores output by IUPred
                            TotalMeanISD, ISD_RawScore, PfamISDAverages = CalculateIUPred(ProteinSequence, PfamIDs,PfamStarts, PfamStops)
                            # Then, each list that will be uploaded to the MySQL database is converted to a string
                            # so that it can be uploaded successfully
                            PfamIDs = "%s" %PfamIDs
                            PfamStarts = "%s" %PfamStarts
                            PfamStops = "%s" %PfamStops
                            InterProScanEval = "%s" %InterProScanEval
                            ISD_RawScore = "%s" %ISD_RawScore
                            PfamISDAverages = "%s" %PfamISDAverages
                        # If the sequence does not pass the quality control filter, the sequences are assigned
                        # NULL for all three ISD values
                        else:
                            TotalMeanISD = 'NULL'
                            ISD_RawScore = 'NULL'
                            PfamISDAverages = 'NULL'
                            PfamIDs = InterProScanResultsDictionary[CodingID]['PfamUID']
                            PfamStarts = InterProScanResultsDictionary[CodingID]['PfamStart']
                            PfamStops = InterProScanResultsDictionary[CodingID]['PfamStop']
                            InterProScanEval = InterProScanResultsDictionary[CodingID]['Eval']
                            PfamIDs = "%s" %PfamIDs
                            PfamStarts = "%s" %PfamStarts
                            PfamStops = "%s" %PfamStops
                            InterProScanEval = "%s" %InterProScanEval
                        # The species is the name of the speciesDirectory where the file was pulled from
                        species = speciesDirectory
                        repository = 'NCBI'
                        # The UID that was used for each key is then decomposed into its three parts
                        # for upload into the MySQL database
                        GeneID, ProteinID,SequenceAccession = CodingID.split('||')


                        '''
                        The final block of code is used to upload all the results into a MySQL database.

                        The personal information of the user used to connect to the database is stored in a function in the DataExtractionSubmodules file so their information can be kept private. 

                        Each sequence is uploaded immediately following its run through IUPred
                        '''
                        # The information needed to connect to the database is loaded from the submodules file
                        Database, User, Host, Password = LoadConnectionInformation()
                        # The information is then used to connect to the database
                        cnx = mysql.connector.connect(user = User,
                                                      password = Password,
                                                      host = Host,
                                                      database = Database)

                        # a cursor is declared
                        mycursor = cnx.cursor(buffered = True)
                        # An insertion statement is then defined for the coding sequences table
                        sqlUploadCodingDataStatement = "INSERT INTO TestUploadTable_Coding (Repository,SpeciesUID,NewickSpecies,NCBIAccessionID,GeneID,ProteinID,CodingSequence) VALUES (%s,%s,%s,%s,%s,%s,%s)"
                        
                        # The coding upload data is then defined as a tuple
                        uploadCodingData = (repository,SpeciesUID,species,SequenceAccession,GeneID, ProteinID, CodingSequence)
                        
                        # and the commend to upload the data is executed
                        mycursor.execute(sqlUploadCodingDataStatement,uploadCodingData)

                        # The UID auto-generated by the SQL table is then extracted so the protein data can be
                        # associated with the coding genes datatable with a unique ID. The coding sequences and protein sequences
                        # are stored in different datatables because there exists a maximum limit on the length of a row in a MySQL
                        # datatable, and some coding sequences can use the entirety of the data limit
                        extractCodingUIDStatement= "SELECT UID FROM TestUploadTable_Coding WHERE UID = (SELECT MAX(UID) FROM TestUploadTable_Coding)"
                        # The UID for the coding sequence is extracted from where it was just uploaded
                        mycursor.execute(extractCodingUIDStatement)
                        for output in mycursor:
                            codingSeqUID = output[0]
                        # An insertion statement is then defined for uploading data into the protein table
                        sqlUploadProteinDataStatement = "INSERT INTO TestUploadTable_Protein (CodingSeqTableUID, SpeciesUID, NewickSpecies, GeneID, ProteinID, ProteinSequence, PfamUID, PfamStart, PfamStop, InterProScanEval, NoCysMeanISD, NoCysISD_RawIUPredScores, PfamDomainMeanISDValues) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                        # The data that will be uploaded to the protein datatable is then defined
                        uploadProteinData = (codingSeqUID, SpeciesUID,species, GeneID,ProteinID,ProteinSequence, PfamIDs, PfamStarts, PfamStops, InterProScanEval, TotalMeanISD, ISD_RawScore, PfamISDAverages)
                        # And the insertion statement is executed with with the data
                        mycursor.execute(sqlUploadProteinDataStatement, uploadProteinData)
                        # All data is committed to the datatable 
                        cnx.commit()
                        # and the connection is closed
                        cnx.close()

                                                     
# Once the entire program is run, the user is notified of how many species wound up in their final dataset
print('Number of species processed: %s' %SpeciesProcessed)
