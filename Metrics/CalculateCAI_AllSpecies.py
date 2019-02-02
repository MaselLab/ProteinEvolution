import os, sys, json, csv, mysql.connector, datetime

'''
Author : Sara Willis
Date   : February 1, 2019
--------------------------

The purpose of this file is to determine the codon adaptation indices (CAI) for individual species. This is a metric that describes codon usage patterns. 

CAI is calculated in the following way:

   1) First an RSCU value has to be calculated. There is one RSCU value for each of the 64 codons. It is defined as the following: 

           X_i/[ 1/n * Sum from 1 to n of X_j] 

                        n
                      -----        X_i
                      \        -----------
           RSCU_i  =   \            N
                       /           ---
                      /         1  \  
                      -----     -  /     X_j
                      i = 1     N  ---
                                  j = 1

      where X_i is the codon in question, and all X_j's are the codons that code for that same amino acid. n is the degeneracy of that codon, i.e. the number of codons that code for that specific amino acid. 

   2) Once the RSCU values have been calculated, the relative adaptedness, or w_i, needs to be found. w_i is defined as:

           RSCU_i/RSCU_max
  
      RSCU_i is the RSCU value for codon i, RSCU_max is the maximum RSCU value associated with all codons synonymous to codon i

   3) The CAI is then defined as the geometric mean of the relative adaptedness values:

           CAI = (product from 1 to L of w_i)^1/L

               ---            ---  1/L 
               |      L         | 
               |    -----       |
               |     | |        |
         CAI = |     | |   w_i  |
               |    i = 1       |
               ---            ---

      Where L is the number of codons in the protein-coding sequences 


The sequences are extracted from a SQL database one species at a time, are then used to calculate the CAI value for that species, and the results are uploaded into the database
'''

####################################################################################
#                                 User Information                                 #
####################################################################################

# SQL Connection Information
Database = ''
User = ''
Host = ''
Password = ''

# The tables where the coding sequences are stored
# ** These are stored as a list so that we can iterate through them **
CodingTables = ['EnsemblGenomes_Coding_Complete', 'NCBIGenomes_Coding_Complete']



####################################################################################
#                             Program Executes Below                               #
####################################################################################

# A connection to the SQL server is established
cnx = mysql.connector.connect(user = User,
                              password = Password,
                              host = Host,
                              database = Database)
mycursor = cnx.cursor(buffered = True)

# We'll pull species from all coding data tables one at a time
for CodingTable in CodingTables:
    # We find the min/max value of the species UIDs in the table and then
    # iterate through all integers in that range so we only deal with one
    # species at a time
    SelectMaxSpeciesUID = "SELECT MAX(SpeciesUID) FROM %s"%CodingTable
    mycursor.execute(SelectMaxSpeciesUID)
    MaxSpeciesUID = mycursor.fetchone()[0]

    SelectMinSpeciesUID = "SELECT MIN(SpeciesUID) FROM %s"%CodingTable
    mycursor.execute(SelectMinSpeciesUID)
    MinSpeciesUID = mycursor.fetchone()[0]
    
    for i in range(MinSpeciesUID,MaxSpeciesUID):
        print('Extracting SpeciesUID: %s'%i)
        totalCodonCount = 0
        # We'll keep dictionaries of all the codons that we count, their RSCU values,
        # and their relative adaptedness values, sorted by which amino acid they
        # correspond to. This will make calculating the CAI relatively easy and painless
        RawCount = {'F':{'TTT':0,'TTC':0},
                    'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                    'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0},
                    'Y':{'TAT':0,'TAC':0},
                    '*':{'TAA':0,'TAG':0,'TGA':0},
                    'C':{'TGT':0,'TGC':0},
                    'W':{'TGG':0},
                    'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                    'H':{'CAT':0,'CAC':0},
                    'Q':{'CAA':0,'CAG':0},
                    'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0},
                    'I':{'ATT':0,'ATC':0,'ATA':0},
                    'M':{'ATG':0},
                    'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                    'N':{'AAT':0,'AAC':0},
                    'K':{'AAA':0,'AAG':0},
                    'S':{'AGT':0,'AGC':0},
                    'R':{'AGA':0,'AGG':0},
                    'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                    'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                    'D':{'GAT':0,'GAC':0},
                    'E':{'GAA':0,'GAG':0},
                    'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

        RSCUTable = {'F':{'TTT':0,'TTC':0},
                     'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                     'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0},
                     'Y':{'TAT':0,'TAC':0},
                     '*':{'TAA':0,'TAG':0,'TGA':0},
                     'C':{'TGT':0,'TGC':0},
                     'W':{'TGG':0},
                     'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                     'H':{'CAT':0,'CAC':0},
                     'Q':{'CAA':0,'CAG':0},
                     'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0},
                     'I':{'ATT':0,'ATC':0,'ATA':0},
                     'M':{'ATG':0},
                     'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                     'N':{'AAT':0,'AAC':0},
                     'K':{'AAA':0,'AAG':0},
                     'S':{'AGT':0,'AGC':0},
                     'R':{'AGA':0,'AGG':0},
                     'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                     'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                     'D':{'GAT':0,'GAC':0},
                     'E':{'GAA':0,'GAG':0},
                     'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

        RelativeAdaptednessTable = {'F':{'TTT':0,'TTC':0},
                                    'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
                                    'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0},
                                    'Y':{'TAT':0,'TAC':0},
                                    '*':{'TAA':0,'TAG':0,'TGA':0},
                                    'C':{'TGT':0,'TGC':0},
                                    'W':{'TGG':0},
                                    'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
                                    'H':{'CAT':0,'CAC':0},
                                    'Q':{'CAA':0,'CAG':0},
                                    'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0},
                                    'I':{'ATT':0,'ATC':0,'ATA':0},
                                    'M':{'ATG':0},
                                    'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
                                    'N':{'AAT':0,'AAC':0},
                                    'K':{'AAA':0,'AAG':0},
                                    'S':{'AGT':0,'AGC':0},
                                    'R':{'AGA':0,'AGG':0},
                                    'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
                                    'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
                                    'D':{'GAT':0,'GAC':0},
                                    'E':{'GAA':0,'GAG':0},
                                    'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}}

        # We start by attempting to extract a set of coding sequences sharing the specified species UID
        SelectionStatement = "SELECT CodingSequence FROM %s WHERE SpeciesUID = %s"%(CodingTable,i)
        mycursor.execute(SelectionStatement)
        speciesSpecificResults = mycursor.fetchall()
        # If the species UID doesn't exist in the table, we move on, otherwise we iterate through all
        # extracted coding sequences to count the occurences of each codon
        if speciesSpecificResults != []:
            for CodingSequence in speciesSpecificResults:

                CodingSequence = CodingSequence[0]
                sequenceLength = len(CodingSequence)
                # We partition each coding sequence into a list codons and then compare them to our dictionary
                codonList = [CodingSequence[n:n+3] for n in range(0,sequenceLength,3)]
                # There are some "codons" that we will be ignoring. Specifically, anything that has an x or some other
                # irregular character in it
                for AA in RawCount:
                    for Codon in RawCount[AA]:
                        CodonCount = codonList.count(Codon)
                        RawCount[AA][Codon] += CodonCount
                        totalCodonCount += CodonCount
            # Once the raw counts are found for each codon, we then start calculating RSCU_i values (see the beginning of
            # this script for the specifics). 
            for AA in RawCount:
                Sum = sum(RawCount[AA].values())
                for Codon in RawCount[AA]:
                    # In some weird cases an amino acid never shows up in a particular genome... if this is the case, we set
                    # the RSCU_i value to one for all codons corresponding to that amino acid. We do this because later we'll
                    # take the geometric mean of all relative adaptedness values and multiplying by one will not affect this value
                    # In essence, setting them to one "silences" these.
                    if Sum != 0:
                        RSCU = RawCount[AA][Codon]/Sum
                    else:
                        RSCU = 1
                    RSCUTable[AA][Codon]=RSCU
            # Once all RSCU_i values are found, we can calculate the relative adaptedness values and save those to our dictionary
            for AA in RSCUTable:
                MaxRSCU = max(RSCUTable[AA].values())
                for Codon in RSCUTable[AA]:
                    RelativeAdaptedness = RSCUTable[AA][Codon]/MaxRSCU
                    RelativeAdaptednessTable[AA][Codon] = RelativeAdaptedness

            # Finally, we are able to calculate the CAI value
            
            # We start by setting the CAI to one
            CAI = 1
            # We then go through the raw codon count dictionary and use those values to calculate the CAI (see beginning of script for
            # details)
            for AA in RawCount:
                for Codon in RawCount[AA]:
                    CodonCount = codonList.count(Codon)
                    if CodonCount != 0:
                        # We raise to the power of CodonCount/totalCodonCount rather than (relativeAdaptedness^CodonCount)^(1/L)
                        # because the latter sometimes forces the value to become so small that it exceeds python's capabilities
                        # and becomes 1 before it raises it to the power of 1/L, which forces the entire CAI to become 0
                        CAI = CAI * ((RelativeAdaptednessTable[AA][Codon])**(CodonCount/totalCodonCount))
            # Once the CAI is found, the result is uploaded to the table in the SQL database
            UpdateStatement = "UPDATE SpeciesList SET CAI=%s WHERE SpeciesUID = %s"%(CAI,i)
            mycursor.execute(UpdateStatement)
            cnx.commit()
            print('Updated!\n\n')
cnx.close()
