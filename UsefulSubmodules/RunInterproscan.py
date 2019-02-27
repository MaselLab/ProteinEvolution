'''

Author : Sara Willis
Date   : Wednesday February 27, 2019


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
import os, csv

def RunInterProScan(FastaFilename, InterProScanExecutablePath, InterProScanOutputPath):
    # InterProScan is run only using Pfam as the reference database using the input specifications.
    # Output is redirected to a tab-delimited file.
    os.system('%s --applications Pfam --input %s --outfile %s -f TSV' %(InterProScanExecutablePath,FastaFilename,InterProScanOutputPath))

    # The results dictionary is defined
    InterProScanResultsDictionary = {}
    with open(InterProScanOutputPath, 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            ProteinAccession = row[0]
            SequenceMD5Digest = row[1]
            SequenceLength = row[2]
            Analysis = row[3]
            SignatureAccession = row[4]
            SignatureDescription = row[5]
            StartLocation = row[6]
            StopLocation = row[7]
            Score = row[8]
            Status = row[9]
            Date = row[10]

            if ProteinAccession not in InterProScanResultsDictionary:
                InterProScanResultsDictionary[ProteinAccession] = {'PfamUID':[SignatureAccession], 'PfamStart':[StartLocation], 'PfamStop':[StopLocation], 'Eval': [Score]}
            else:
                InterProScanResultsDictionary[ProteinAccession]['PfamUID'].append(SignatureAccession)
                InterProScanResultsDictionary[ProteinAccession]['PfamStart'].append(StartLocation)
                InterProScanResultsDictionary[ProteinAccession]['PfamStop'].append(StopLocation)
                InterProScanResultsDictionary[ProteinAccession]['Eval'].append(Score)
    if os.path.exists(InterProScanOutputPath) == True:
        os.remove(InterProScanOutputPath)
    return InterProScanResultsDictionary




