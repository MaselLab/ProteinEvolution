import sys, json, os, csv, time, shutil, ftputil, string, re, datetime, mysql.connector
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein



#################################################################################
########################## LoadConnectionInformation ############################
#################################################################################

'''
Because information will be uploaded into a personal MySQL database, if the user wants to open the main script to edit, they may wish to keep their login information private. This function will load usernames and passwords to keep them hidden.
'''

def LoadConnectionInformation():
    # Name of MySQL database
    Database = ''
    # Username
    user = ''
    # Database IP
    host = ''
    # User's MySQL password
    password = ''
    
    return Database, user, host, password


#################################################################################




#################################################################################
################################ CreateFastaFile ################################
#################################################################################

'''
The purpose of this function is to create a fasta file with sequences designated by the user. The function taskes in three arguments:

1) UID will be the identifier for the sequence in the fasta file
2) Protein is the protein sequence
3) HitsFastaFile is the filename the user wishes the file to be saved as

IMPORTANT:

The function writes to the file using the "append" function, so sequences can be written to the file by calling the function in a for loop. This, however, means that the user must remember to delete old files as, if they are not, new data may be inadvertently appended to the ends of old files
'''

def CreateFastaFile(UID, Protein,HitsFastaFile, Description='UnknownDescription'):
    # we save sys.stdout as a new variable so that we can redirect what's printed
    # to the terminal to a file, and then revert to printing to the terminal at the
    # end of the function's use
    orig_stdout = sys.stdout
    # BioPython is used to generate a sequence record so that it's module SeqIO
    # can be called to write in fasta format
    record = SeqRecord(Seq(Protein), id= '%s'%UID, description = Description)
    # stdout is redirected to a file with the name chosen by the user
    sys.stdout = open('%s'%HitsFastaFile, 'a')
    print(record.format('fasta'))
    # Once it has been printed, sys.stdout is closed and the program will start
    # printing to the terminal again
    sys.stdout.close()
    # Future output is then redirected to the terminal
    sys.stdout=orig_stdout

#################################################################################



#################################################################################
################################ RemoveAllFiles  ################################
#################################################################################

'''
This function is a simple way of deleting all old temporary files. This is done to ensure that in an iterative program, no old 
'''

def RemoveAllFiles(TempFilesDirectoryPath):
    # Checks to see if TempFiles directory exists
    if os.path.exists(TempFilesDirectoryPath) == True:
        # if it does, it is removed along with all its contents
        shutil.rmtree(TempFilesDirectoryPath)
    # An empty TempFiles directory is then created
    os.mkdir(TempFilesDirectoryPath)
    # a short pause is then used to ensure all files are deleted
    time.sleep(0.1)

#################################################################################



#################################################################################
################################ Calculate IUPred ###############################
#################################################################################

'''
This function is used to caluclate the intrinsic structural disorder of a protein sequence and sections of the sequence based on Pfam coordinates provided by the user. 

These analyses are done using the program IUPred*. The executable for the program as well as all relevant files should be located in the same directory as this script.

All analyses are done on sequences with all cysteines removed. 

=================================================================
INPUT
-----
     Protein -- string format
     PfamIDs -- The UIDs for each Pfam associated with the provided gene. There can be multiple Pfams associated with each gene and should be input as a list
     PfamStart/Stop -- The starting/ending location of each Pfam domain in protein-space. Should be in list format. The indices of each start/stop location should correspond with the index of the Pfam UID it is associated with.

Output
------
     totalMeanISD -- The average of all raw scores produced from IUPred for the entire cysteineless protein sequence. Output as a float
     ISD_RawScore -- A list with all the raw scores output from IUPred
     PfamISDAverages -- A list of average ISD values calculated for each Pfam domain. The indices of the averages in the output list will correspond to the indices of the PfamUIDs they correspond with. They are caluclated by using the raw scores and averaging over the regions associated with each Pfam domain.

=================================================================

*IUPred: web server for the prediction of intrinsically unstructured regions of proteins based on estimated energy content
Zsuzsanna Dosztányi, Veronika Csizmók, Péter Tompa and István Simon 
Bioinformatics (2005) 21, 3433-3434. 
'''
def CalculateIUPred(ProteinSequence, PfamIDs, PfamStarts, PfamStops):
    # Defines a filename for a temporary fasta file for IUPred to use
    TempFastaFile = 'Temp.fasta'
    # A temporary file for the output from IUPred to be printed to 
    IUPredOutputFile = 'IUPredOutputFile.txt'
    # A list of characters used to parse the space-delimited IUPred output file
    unwanted = ['',"'",' ']
    # A new sequence is made by removing all cysteines from the protein sequence
    NoCysProteinSequence = ProteinSequence.replace('C','')
    # Creates a temporary Fasta file for Iupred to read. If a file with that name
    # already exists, it is removed to avoid printing multiple sequences to a file
    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)
    # The protein sequence with no cysteines included is printed to a
    # temporary fasta file
    CreateFastaFile('Unknown',NoCysProteinSequence,TempFastaFile)
    # If an IUPred output file already exists, it's removed
    if os.path.exists(IUPredOutputFile) == True:
        os.remove(IUPredOutputFile)
        
    # Executes IUPred and prints the output to a temporary file
    os.system('./iupred %s long >%s' %(TempFastaFile, IUPredOutputFile))
    # A string where the raw ISD scores will be stored is defined
    ISD_RawScore = []
    # The temporary output IUPred file is opened to be read
    with open(IUPredOutputFile,'r') as f:
        # Each line is read as a list
        reader = csv.reader(f)
        for row in reader:
            # The line is redefined as a string so it can be split back into
            # list using spaces as a delimited
            splitRow = row[0].split(' ')
            # a new row is defined only including characters that are relevant
            # (not included in the list 'unwanted'
            newRow = []
            for element in splitRow:
                if element not in unwanted:
                    newRow.append(element)
            # If a new run starts with the character '#', it is discarded
            # (these correspond to human-readable information and do not
            # contain relevant data)
            if '#' in newRow:
                pass
            # Otherwise, there are three relevant values to extract
            else:
                # 1) The amino acid number in the full protein
                AANumber = int(newRow[0])
                # 2) The amino acid in the sequence
                AA = newRow[1]
                # 3) The calculated ISD value (0 <= ISDValue <= 1)
                ISDValue = float(newRow[2])
                # All ISD values are appended to a list to return to the user
                ISD_RawScore.append(ISDValue)
    # A list of average ISD values for each Pfam domain is defined
    PfamISDAverages = []
    # Then each Pfam is used to find the ISD values associated with it
    for PfamID in PfamIDs:
        # The index of the Pfam is determined to determine its
        # starting/stopping locations
        PfamIndex = PfamIDs.index(PfamID)
        PfamStart = int(PfamStarts[PfamIndex])
        PfamStop = int(PfamStops[PfamIndex])
        # Because the start/stop locations are associated with the protein
        # including its cysteines the indices need to be transformed into
        # cysteineless protein space. First, the Pfam peptide sequence is
        # determined using the start/stop locations
        PfamPeptide = ProteinSequence[PfamStart:PfamStop]
        # The cysteines are then removed from the Pfam peptide sequence
        NoCysPfamPeptide = PfamPeptide.replace('C','')
        # The starting index of this new subsequence within the full
        # cysteineless sequence is determined
        newPfamStart = NoCysProteinSequence.index(NoCysPfamPeptide)
        # and the end is determined by adding the new start with the
        # length of the subsequence
        newPfamStop = newPfamStart + len(NoCysPfamPeptide)
        # All raw scores associated with this subsequence are then pulled
        # from the list of raw IUPred scores
        PfamISDList = ISD_RawScore[newPfamStart:newPfamStop]
        # and are averaged to give the mean ISD for that protein's Pfam
        # domain
        PfamISDAverages.append(np.mean(PfamISDList))
        
        
    # The total mean ISD is the mean of all the raw scores for the whole protein 
    totalMeanISD = np.mean(ISD_RawScore)
    # The temporary files are then removed
    if os.path.exists(TempFastaFile) == True:
        os.remove(TempFastaFile)
    if os.path.exists(IUPredOutputFile) == True:
        os.remove(IUPredOutputFile)
    # The user is notified that IUPred ran successfully and the time of completion
    # and the relevant information is returned to the user.
    return totalMeanISD, ISD_RawScore, PfamISDAverages
        

#################################################################################


#################################################################################
################################ RunInterProScan ################################
#################################################################################

'''
This function will take in a fasta file containing protein sequences and will run them through InterProScan to determine the PfamIDs they are associated with as well as the expectation values

===============================================

INPUT
-----
     1) Fasta file containing the protein sequences to search against the database
     2) Path to InterProScan executable (interproscan.sh). This should be relative to the directory where this script is located
     3) Output Filename and path. This should be relative to the InterProScan directory, and should be directed to the 

OUTPUT
------
     Nested dictionary with the following format:
     {ProteinAccession: {PfamUID: [PfamUID1, ..., PfamUIDn], PfamStart: [PfamUID1_start, ..., PfamUIDn_start], PfamStop: [PfamUID1_stop, ..., PfamUIDn_stop], Eval: [PfamUID1_eval, ..., PfamUIDn_eval]}

===============================================

'''

def RunInterProScan(FastaFilename, InterProScanExecutablePath, InterProScanOutputPath):
    # InterProScan is run only using Pfam as the reference database using the input specifications.
    # Output is redirected to a tab-delimited file.
    os.system('%s --applications Pfam --input %s --outfile %s -f TSV' %(InterProScanExecutablePath,FastaFilename,InterProScanOutputPath))

    # The results dictionary is defined
    InterProScanResultsDictionary = {}
    # The InterProScan output file is opened for parsing
    with open(InterProScanOutputPath, 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        # Each row is parsed for the relevant data
        for row in reader:
            # Accession of the sequence that was searched against the Pfam database
            ProteinAccession = row[0]
            # Information about the protein. Will not be used
            SequenceMD5Digest = row[1]
            # Length of the searched sequence
            SequenceLength = row[2]
            # More descriptions that will not be used
            Analysis = row[3]
            # Pfam hit ID
            SignatureAccession = row[4]
            # Description of Pfam
            SignatureDescription = row[5]
            # Location of Pfam domain start in protein coordinates
            StartLocation = row[6]
            # Location of Pfam domain stop in protein coordinates
            StopLocation = row[7]
            # Evalue associated with Pfam hit
            Score = row[8]
            # More descriptions that won't be used
            Status = row[9]
            # Date -- will not be used
            Date = row[10]

            '''
            Each entry in the InterProScanResultsDictionary will have nested dictionary with the following format:

            {Protein Accession: {PfamUID: [Pfam_1, ..., Pfam_n], PfamStart: [Pfam_1_start, ..., Pfam_n_start], PfamStop: [Pfam_1_stop, ..., Pfam_n_stop], Eval: [Pfam_1_eval, ..., Pfam_n_eval]}}
            '''
            # If the searched protein is not already in the dictionary, it is added along with all
            # the relevant information 
            if ProteinAccession not in InterProScanResultsDictionary:
                # Each entry has its relevant information added as a length-one list so that further Pfam
                # data can be appended to the lists
                InterProScanResultsDictionary[ProteinAccession] = {'PfamUID':[SignatureAccession], 'PfamStart':[StartLocation], 'PfamStop':[StopLocation], 'Eval': [Score]}
            else:
                # If the searched protein is already in the dictionary, then the relevant information is
                # appended to the preexisting lists
                InterProScanResultsDictionary[ProteinAccession]['PfamUID'].append(SignatureAccession)
                InterProScanResultsDictionary[ProteinAccession]['PfamStart'].append(StartLocation)
                InterProScanResultsDictionary[ProteinAccession]['PfamStop'].append(StopLocation)
                InterProScanResultsDictionary[ProteinAccession]['Eval'].append(Score)
    # Once the file is completely parsed, it is removed
    if os.path.exists(InterProScanOutputPath) == True:
        os.remove(InterProScanOutputPath)
    # The user is then returned the dictionary
    return InterProScanResultsDictionary



#################################################################################


#################################################################################
########################### CreateNCBISpeciesList ###############################
#################################################################################

# This function takes in a TSV file that contains the species that need to be
# pulled from the NCBI databases. The species are extracted from the file and
# added to a list that is returned to the user

# The function takes in the filename as input
def CreateNCBISpeciesList(NCBISpeciesFilename):
    # An empty list is then defined
    speciesList = []
    # The file is opened for parsing
    with open(NCBISpeciesFilename, 'r') as f:
        # The reader is defined to read the file as a TSV
        reader = csv.reader(f, delimiter = ',')
        # Each row is then parsed
        for row in reader:
            # The second entry in the file contains the species names
            species = row[1]
            # Each species name is added to the defined list
            speciesList.append(species)
    # Once the file has been completely parsed, the list is returned to the user
    return speciesList

#################################################################################


#################################################################################
################################ DownloadFTPFile  ###############################
#################################################################################

# The purpose of this function is to download a declared file from NCBI from their
# FTP service
def Download(SourceFilename):
    # The host is declared
    host = ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')
    # The genbank directory is navigated to
    host.chdir('/genbank')
    # The file is then downloaded with the same name as appears in the directory
    host.download(SourceFilename,SourceFilename)
    # The zipped file is then unzipped
    os.system('gzip -d %s' %SourceFilename)
    # And the function exits
    return

#################################################################################



#################################################################################
######################## Coding Sequence Quality Control  #######################
#################################################################################
'''
The purpose of this function is to check a coding sequence to make sure that it:

    1) Starts with a start codon
    2) Doesn't contain any in-frame stop codons
    3) Is a multiple of 3 nucleotides in length.

The function takes a coding sequence as input and returns either:

    1) True if the coding sequence passes all filters
    or 
    2) False if any of the conditions are not met
'''

def CodingSequenceQualityControl(CodingSequence):
    # First, the coding sequence needs to start with a start codon
    if CodingSequence[:3] != 'ATG':
        # if it doesn't, the user is notified the sequence did not pass the quality filters
        return False
    # If the coding sequence isn't a multiple of 3 in length, it fails the quality check
    elif len(CodingSequence)%3 != 0:
        return False
    else:
        # Lastly, in-frame stop codons are searched for
        # The last three nucleotides are excluded since we don't want to throw away a
        # sequence for ending with a stop codon
        for i in range(0, len(CodingSequence)-3, 3):
            Codon = CodingSequence[i:i+3]
            # If any stop codons are found in-frame in the body of the sequence, the
            # quality check is failed
            if Codon == 'TAA' or Codon == 'TGA' or Codon =='TAG':
                return False
            else:
                # If the sequence passes all filters, it passes
                return True

#################################################################################


#################################################################################
############################## Attempt File Download  ###########################
#################################################################################

def AttemptFileDownload(directory, subdirectory, genomeType, TempFilesDirectoryPath):
    host = ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')
    # Change to the indicated genomeType
    host.chdir('/genomes/refseq/%s/%s/%s'%(directory, subdirectory, genomeType))
    print('\nDirectory Found')
    # Determinet the name of the one included subdirectory
    contentOfSpeciesDirectory = host.listdir(host.curdir)[0]
    # navigate into that subdirectory
    host.chdir('/genomes/refseq/%s/%s/%s/%s'%(directory, subdirectory,genomeType,contentOfSpeciesDirectory))
    # determine if the directory has coding sequences associated with it
    filesList = host.listdir(host.curdir)
    CodingSequenceFileExists = False
    for File in filesList:
        # Checks each file/directory in the current location for a file
        # with the extension that indicates the coding sequence file
        if '_cds_from_genomic.fna.gz' in File:
            CodingSequenceFile = File
            # If the file found, the loop is exited
            CodingSequenceFileExists = True
            break
        # If the coding sequence file exists, it is downloaded
    if CodingSequenceFileExists == True:
        print('\nDownloading File...\n')
            
        # The file is saved in the temporary files directory
        TargetFilename = TempFilesDirectoryPath + CodingSequenceFile
        # Downloads the file
        host.download(CodingSequenceFile, TargetFilename)
        # Checks to see that the file was downloaded
        if os.path.exists(TargetFilename) == True:
            print('\nFile Successfully Downloaded.\n')
            print('\nUnzipping File...\n')
            # if the file was downloaded correctly, it is unzipped
            os.system('gzip -d %s' %TargetFilename)
            # Checks to see that the file was unzipped successfully
            UnzippedTargetFilename = TargetFilename.replace('.gz','')
            # Checks to make sure the file was unzipped successfully
            if os.path.exists(UnzippedTargetFilename) == True:
                print('\nFile Successfully Unzipped.\n')
        else:
            # If the coding sequence file doesn't exist, the function prints
            # that there was an issue and shuts down
            print('There was a problem downloading the file')
        return UnzippedTargetFilename
                    

#################################################################################


#################################################################################
############################## Attempt File Download  ###########################
#################################################################################

'''
The purpose of this function is to determine whether or not a species is present in the primary species table on MySQL. If it's not, it is uploaded to the table and the Sequence UID that is associated with it is extracted and returned to the user. If it already is in the table, then that Sequence UID is extracted for the user and returned.
'''

# The function takes the species name in as input
def FormatAndExtractSpeciesUID(species):
    # The connection information is loaded for the user's MySQL database
    Database, User, Host, Password = LoadConnectionInformation()
    # Connects to the database with the loaded information
    cnx = mysql.connector.connect(user = User,
                                  password = Password,
                                  host = Host,
                                  database = Database)
    # A cursor is declared so the database can be interacted with
    mycursor = cnx.cursor(buffered = True)
    # An extraction statement from the temporary species database is defined
    extractSpeciesData = "SELECT * FROM Species_NotInEnsembl_CompleteInGold_InTimetree_NotInSpeciesList"
    # The entire row is extracted for each species
    mycursor.execute(extractSpeciesData)
    # The cursor is then parsed
    for speciesData in mycursor:
        # if the input species is found in an entry, the entry is parsed
        if speciesData[1] == species:
            # The original species UID from the temporary database
            oldSpeciesUID = speciesData[0]
            # The species name
            Name = speciesData[1]
            # The status of the study ID in the GOLD database
            Status = speciesData[2]
            # The study ID
            StudyID = speciesData[3]
            # the loop is then exited
            break
    # A statement to select from the permanent species table is defined
    searchSpeciesTable = "SELECT * FROM SpeciesList"
    # The search statement is then executed
    mycursor.execute(searchSpeciesTable)
    # By default, the variable speciesFound is set to False and will only be set
    # to True if the input species is found in the permanent species table
    speciesFound = False
    # The data from the permanent species table is parsed
    for speciesData in mycursor:
        # If the species is found in the permanent table, speciesFound is set to True
        # and the loop is exited
        if speciesData[2] == species:
            # The speciesUID is returned to the user if it's found in the permanent table
            speciesUID = speciesData[0]
            # speciesFound is set to True so species will not be uploaded to permanent table
            speciesFound = True
            # loop is exited
            break
    # If the species wasn't found in the permanent species table, then it will be added
    if speciesFound == False:
        # An insertion statement is declared so that the species can be uploaded to the permanent
        # species datatable
        uploadSpeciesData = "INSERT INTO SpeciesList (NewickSpeciesName, HighestGoldStatus,StudyID, EnsemblAccession,EnsemblDB,InBiomart) VALUES (%s,%s,%s,%s,%s,%s)"
        # The data that need to be uploaded come from the temporary table or are declared, and are
        # put into a tuple to be uploaded
        SpeciesData = (Name, Status,StudyID, 'N/A','NCBI','N/A')
        # The upload statement is then executed with the relevant data
        mycursor.execute(uploadSpeciesData, SpeciesData)
        # The species UID is then extracted from where the species was just uploaded for the user
        extractSpeciesUID = "SELECT SpeciesUID FROM SpeciesList WHERE SpeciesUID = (SELECT MAX(SpeciesUID) FROM SpeciesList)"
        mycursor.execute(extractSpeciesUID)
        # The Species UID is then pulled from the cursor (which is returned as a length one tuple)
        for uid in mycursor:
            speciesUID = uid[0]
    # all changes are committed to the database
    cnx.commit()
    # and the connection is closed
    cnx.close()
    # And the user is returned the speciesUID extracted from the permanent datatable
    return speciesUID
                                  


#################################################################################
