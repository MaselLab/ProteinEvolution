from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys, json, csv, re, string, mysql.connector, datetime, random
import numpy as np

'''
Author : Sara Willis
Date   : March 8, 2019

This script is designed to read align pfams all sharing a common ID from many different species. It goes about this by accessing protein data tables (in our dataset there's one from ensembl and one from NCBI) and extracting protein sequences, the pfams that are contained within them protein sequence, the indices where they occur, and a quality check. Proteins that pass the quality check move forward and have the peptide sequences associated with their pfams extracted. 

There is the poosibility that a pfam will show up more than once in a given protein. If this is the case, the scripts looks at the lengths of the pfams and selects one at random from the upper 50th percentile of given lengths.

Once all the pfams meeting the previously described conditions, the script goes about the process of aligning all peptides associated with a particular pfam UID using Clustal Omega and uploads each aligned sequence to a MySQL data table. 

In some cases, Clustal Omega is unable to align all sequences associated with a particular pfam because there are too many. Thus far, it's not clear what to do with those. Another alignment program may be necessary. Until that's implemented, the script keeps track of these problematic pfams and prints them to a flat file for the user to look at once the script has been run.

DEPENDENCIES
------------

This script requires Clustal Omega to run which can be downloaded as a precompiled executable from: 

http://www.clustal.org/omega/#Download

The script also requires the non-standard python module mysql.connector
'''
########################################################################
#                        User-Specific Information                     #                     
########################################################################

# User's MySQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# Databases where the sequences are stored. We need to access whether the "host" protein
# has passed quality filters, so we include the protein metrics table (where that information
# is stored)
MulticellularProteinTable = "Genomes_Multicellular_Protein"
MulticellularProteinMetricsTable = "Genomes_Multicellular_ProteinMetrics"
NCBIProteinTable = "NCBIGenomes_Protein_Complete"
NCBIProteinMetricsTable = "NCBIGenomes_ProteinMetrics_Complete"

PfamAlignmentsTable = 'PfamAlignments_copy_copy'

# Clustal is greedy and will use every available core it can find. You can teach it
# restraint 
numberOfClustalThreads = 20

# How "chatty" you'd like the script to be
Verbose = True

########################################################################
#                              Submodules                              #                     
########################################################################

# This function creates a temporary fasta file in the current directory
# The user should be careful to delete old files they no longer need when
# using this submodule since it uses the append function instead of overwriting
# preexisting files

def CreateFastaFile(UID, Protein, Description,HitsFastaFile):
    # we save sys.stdout as a new variable so that we can redirect what's printed
    # to the terminal to a file, and then revert to printing to the terminal at the
    # end of the function's use
    orig_stdout = sys.stdout
    # BioPython is used to generate a sequence record so that its module SeqIO
    # can be called to write in fasta format
    record = SeqRecord(Seq(Protein), id= '%s'%UID,description=Description)
    # stdout is redirected to a file with the name chosen by the user
    sys.stdout = open('%s'%HitsFastaFile, 'a')
    print(record.format('fasta'))
    # Once it has been printed, sys.stdout is closed and the program will start
    # printing to the terminal again
    sys.stdout.close()
    sys.stdout=orig_stdout

# This function takes in a dictionary with the following form:
#
# {key : [(value1_1,value2_1), ... , (value1_n,value2_n)]}
#
# and selects a random tuple such that value2_n-value1_n is in the upper 50th
# percentile of the set of all such differences

def ChooseFromUpper50thPercentile(dictionary,key):
    # We first assemble a set of all the differences
    lengths = [k[1]-k[0] for k in dictionary[key]]
    
    if len(lengths) != 1:
        FiftiethPercentile = np.percentile(lengths,50)
        # We then whittle down the set of acceptable values and choose one at random
        AcceptableLengths = [(k[0],k[1]) for k in dictionary[key] if k[1]-k[0] >= FiftiethPercentile]
        randomChoice = random.choice(AcceptableLengths)
    else:
        # There's not much to do if there's only one tuple to deal with, so we just return
        # it to the user
        randomChoice = dictionary[key][0]
    return randomChoice


########################################################################
#                        Program Executes Below                        #                     
########################################################################
    
start_time = datetime.datetime.now()

# Connect to the database
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)



# We want to create a species name dictionary so we can attribute sequences in alignments
# to the species they came from 
SpeciesNameDictionary = {}
mycursor.execute("SELECT NewickSpeciesName,SpeciesUID FROM SpeciesList")
SpeciesResults = mycursor.fetchall()
for SpeciesResult in SpeciesResults:
    NewickSpeciesName, SpeciesUID = SpeciesResult
    SpeciesNameDictionary[SpeciesUID] = NewickSpeciesName

# We need to create a dictionary where we can store and access our pfam sequences
ProteinDictionary = {}

# We need to search through two data tables for each repository so we can grab our filter values
ProteinTables = [(NCBIProteinTable,NCBIProteinMetricsTable),(MulticellularProteinTable,MulticellularProteinMetricsTable)]

for ProteinTable in ProteinTables:
    protein = ProteinTable[0]
    proteinMetrics = ProteinTable[1]
    if Verbose == True:
        print('\n\nExtracting protein sequences and pfam data from %s\n'%ProteinTable[0])
        sys.stdout.flush()
        
    # We use the INNER JOIN function in MySQL so we can grab elements from multiple tables
    ExtractionStatement = "SELECT "+protein+".ProteinSequence,"+protein+".PfamUID,"+protein+".PfamStart,"+protein+".PfamStop,"+protein+".SpeciesUID,"+proteinMetrics+".PassedIUPredFilters FROM "+protein+" INNER JOIN "+proteinMetrics+" ON "+protein+".UID="+proteinMetrics+".ProteinTableUID"
    mycursor.execute(ExtractionStatement)
    ProteinResults = mycursor.fetchall()
    if Verbose == True:
        print('Data Extracted!\nTime Taken: %s\n\nCompiling Pfam Sequence Dictionary'%(datetime.datetime.now()))
        sys.stdout.flush()
        intermediate_time = datetime.datetime.now()


    for result in ProteinResults:
        ProteinSequence,PfamUIDs,PfamStart,PfamStop,SpeciesUID,PassedFilters = result
        if int(PassedFilters) == int(True):
            SpeciesName = SpeciesNameDictionary[SpeciesUID]

            # The pfam UIDs and their starting/stopping indices in the pfam tables are stored as
            # comma-delimited strings. They're converted into lists for easy parsing
            PfamUIDs = [i for i in PfamUIDs.split(',')]
            PfamStart = [int(i) for i in PfamStart.split(',')]
            PfamStop = [int(i) for i in PfamStop.split(',')]

            # Our job now is to collect the sequences from each of the pfams that show up in a given protein
            TempPfamDictionary = {}
            for i in range(0,len(PfamStart)):
                # Pfams can show up in a given protein multiple times. If this happens, we want to cluster
                # the sequences that correspond to that Pfam. We will then select one from the upper 50th percentile
                # of the lengths at random
                if PfamUIDs[i] not in TempPfamDictionary:
                    TempPfamDictionary[PfamUIDs[i]] = [(PfamStart[i],PfamStop[i])]
                else:
                    TempPfamDictionary[PfamUIDs[i]].append((PfamStart[i],PfamStop[i]))

            for PfamUID in TempPfamDictionary:
                # We choose a pfam at random from the upper 50th percentile of the pfam lengths
                start,stop = ChooseFromUpper50thPercentile(TempPfamDictionary,PfamUID)
                # And do a quick quality check
                if stop-start > 0 and stop <= len(ProteinSequence):
                    if PfamUID not in ProteinDictionary:
                        ProteinDictionary[PfamUID] = [(ProteinSequence[start:stop],SpeciesName)]
                    else:
                        ProteinDictionary[PfamUID].append((ProteinSequence[start:stop],SpeciesName))



if Verbose == True:
    print('Pfam Sequence Dictionary Compiled\nTime Taken: %s\n'%(datetime.datetime.now()-intermediate_time))
    sys.stdout.flush()
    intermediate_time = datetime.datetime.now()

# The problematic pfams that are too large to align are saved to a file so the user can check them once
# the script is complete
ProblemPfamsFile = open('PfamsTooLargeToAlign.csv','w')

# Once the protein dictionary has been compiled (and this really should be changed to peptide dictionary
# or pfam dictionary), each pfam has all of its peptide sequences printed to a fasta file for alignment
for pfam in ProteinDictionary:
    if Verbose == True:
        print('Creating fasta file for %s\n'%pfam)
        sys.stdout.flush()
    # If there are more than one sequence, we stick them in a fasta file
    if len(ProteinDictionary[pfam]) != 1:
        for entry in ProteinDictionary[pfam]:
            protein = entry[0]
            species = entry[1]
            CreateFastaFile(pfam,protein,species,'%s.fa'%pfam)
        if Verbose == True:
            print('Aligning sequences\n')
            sys.stdout.flush()
        # Once the fasta file has been created, we attempt to align it using clustal omega
        # force overwrites any previously existing alignment file with the name the program is attempting
        # to write. We also redirect the output to a junk file so it doesn't clutter up our own output file
        os.system('./clustalo --force --threads=%s --infile %s.fa --outfile %s.algn >Clustal.out'%(numberOfClustalThreads,pfam,pfam))
        os.remove('%s.fa'%pfam)
        try:
            # Once the alignment file is created, we extract each aligned sequence and stick it in a MySQL table
            for record in SeqIO.parse('%s.algn'%pfam,'fasta'):
                Sequence = str(record.seq)
                Species = record.description.split(' ')[1]
                PfamUID = record.id
                UploadPfamStatement = "INSERT INTO "+PfamAlignmentsTable+" (PfamUID,NewickSpeciesName,AlignedPeptide) VALUES (%s,%s,%s)"
                UploadPfamValues = tuple((PfamUID,Species,Sequence))
                mycursor.execute(UploadPfamStatement,UploadPfamValues)

            os.remove('%s.algn'%pfam)
            cnx.commit()
        except:
            # If it doesn't exist, that's because Clustal Omega couldn't align the sequences. The program
            # writes the pfam UID to a file for the user to investigate after the program has run
            ProblemPfamsFile.write('%s\n'%pfam)
    
    # If only one pfam exists, there's no point in aligning anything, so we upload the sequence
    else:
        UploadPfamStatement = "INSERT INTO "+PfamAlignmentsTable+" (PfamUID,NewickSpeciesName,AlignedPeptide) VALUES (%s,%s,%s)"
        UploadPfamValues = tuple((pfam,ProteinDictionary[pfam][0][1],ProteinDictionary[pfam][0][0]))
        mycursor.execute(UploadPfamStatement,UploadPfamValues)
        cnx.commit()

cnx.close()
ProblemPfamsFile.close()
if os.path.exists('Clustal.out') == True:
    os.remove('./Clustal.out')
print('Alignments Complete\nTime Taken: %s'%start_time)
