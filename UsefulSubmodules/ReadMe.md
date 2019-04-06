# Useful Submodules

This repository is used to store various useful submodules that can be integrated into other, larger bioinformatics pipelines. 

## R Scripts

### CalculateOptimalLambda.r

This is a short script that is used to determine the optimal lambda value for a box-cox transform 

## Python Scripts


### CodingSequenceQualityControl.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this function is to check a coding sequence to make sure that it:

    1) Starts with a start codon
    2) Doesn't contain any in-frame stop codons
    3) Is a multiple of 3 nucleotides in length.
    
The function takes a coding sequence as input and returns either:

    1) True if the coding sequence passes all filters
    
    or 
    
    2) False if any of the conditions are not met
    
 
### CreateBitScoreMap.py   

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this submodule is to find the heights of the sequence logo of a multiple sequence alignment.
The function accepts a multiple sequence alignment in fasta format as input. Each column in the alignment is used to generate the sequence logo array. This submodule can be used with a change point detection module to find optimal regions to partition a multiple sequence alignment.

For more information on sequence logos and the mathematics used to generate the map: https://en.wikipedia.org/wiki/Sequence_logo
Note: The bit score, when plotted, should look identical to the heights of the sequence logo map that shows up in Geneious for multiple sequence alignments

### CreateFastaFile.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this function is to create a fasta file with sequences designated by the user. The function taskes in three arguments:

 1) The identifier for the sequence in the fasta file
 2) Protein sequence
 3) The saved filename
 
IMPORTANT:

The function writes to the file using the "append" function, so sequence can be written to the file by calling the function in a for loop. This, however, means that the user must remember to delete old files as, if they are not, new data may be inadvertently appended to the ends of old files


### DetectChangePoints.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this function is to read in a list of numerical values and output the detected change points associated with that list. For example, if the user wanted to know how to partition a multiple sequence alignment, they could create a bit score list, feed it into this module, and determine where the alignment should be split. The function has three inputs that need to be specified by the user, as well as one optional parameter:

 [1] The bit scores list is the first argument in the function and is what is analyzed for change points
 [2] penalty is what is used to determine how sensitive the program should be when looking for change points. Too sensitive means every slight variation is counted as a change point, while very strict means some may be missed and only the most extreme changes are counted, so some fiddling may be necessary
 
 [3] minimumSize is the minimum number of elements allowed between change points. This can prevent hundreds of change points showing up in a small bumpy region
 
 [4] The last argument is optional. The user can opt to set this to True if they wish to see a display of the bit scores along with the change points. This can be helpful during testing/debugging so the user can better visualize the data points. 
 
 The output of this function is a list of the coordinates where change points were found. The start and stop locations in the MSA are always counted as change points for convinience. A second list of the averages of the bit scores in each partitioned section of the MSA, divided up by the change points, is also returned to the user. This makes spikes easier to identify. It also allows the user to pick the best change point (most dramatic difference between regions), if the user only wants to focus on one region at a time.
 
Dependencies:

Ruptures: a change point detection module in Python. 
          -- http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html
          -- https://pypi.org/project/ruptures/

### FindAPRs_TangoOutput.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this submodule is to search through a list of raw Tango scores to find runs of aggregation-prone amino acids. The function calculates the density of the AAs in APRs in the sequence to the length of the sequence, as well as the number of AAs in APRs in the sequence. These two values are returned to the user as well as a list of the specific runs. This is in case the user needs to calculate the aggregation scores for subsequences of the protein as well.

***** This differs from Scott's script:
If an unknown amino acid is found within an APR, the sequence is excluded from analyses and the user is returned False. If unknown amino acids exist anywhere else in the protein sequence, so long as they are not in a probable APR, the sequence is retained

**Input:**

      1) AggregationList -- The raw Tango scores in list format. This can be obtained by using the module in the file ReplaceXInTangoOutput.py
      2) AggThreshold -- This is the value a Tango score must exceed to be counted as aggregation prone
      3) RunThreshold -- This is the number of Tango scores exceeding the AggThreshold in sequence for the region to be counted as an APR
      
**Output:**

     A list of lists. Each sublist contains the scores in every aggregation prone region.
     
### KingdomID_function.py

Author: Jenny James
Uploaded: 11th March 2019  

This function groups species in the dataset into their kingdoms, and return a dictionary of kingdoms and species UIDs


### PullHits_HMMEROutput_Tblout.py

Author : Sara Willis
Date   : Wednesday, February 27 2019

This submodule is intended to be used in conjuction with the homology-detection software HMMER. HMMER has the option to output summaries of its runs in tblout format. This submodule can be used to go through the output file to collect hits that match the user-specified evalue cutoffs, plus an additional, optional cutoff that can be used to identify possible false positives. The function will then return to the user a list of all the UIDs associated with the hits passing all filters.
Input:
      1) ECutoff                 -- The max expectation value allowable for a hit to be kept. This Evalue is for the whole sequence
      2) domECutoff              -- an additional Eval cutoff, but for domains only
      3) filename                -- The tblout filename
      4) (Optional) bitBiasRatio -- There are bit scores and bias scores provided by HMMER. HMMER states that sometimes regions of low complexity can
                                    cause low enough Eval scores that a hit is registered as significant. The bias score attempts to correct for this.
                                    The user manual suggests that if the bias score is roughly the same order of magnitude as the bit score, the hit is
                                    possibly a false positive. If the user wishes to use the value bias/bit as a filter, this is supplied as an option


### QualityControlForTango.py

Author : Sara Willis
Date   : February 27, 2019

This function is used to check to see whether a sequence passes imposed quality filters before it's run through the aggregation-prediction software Tango. Specifically, if a protein has three or more unknown residues in sequence, then the protein does not pass the filter and the user is notified to not use the sequence for Tango analysis. 

**Input**: A protein sequence

**Output**: True or False depending on whether the sequence passes the quality check or not.

### ReplaceXInTangoOutput.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this submodule is to take in a file of raw Tango output and return to the user a list of the concatenated raw scores. In the process, it reintroduces unknown amino acids located in the body of the original sequence back into the raw score list as X's

**Input** : 

       1) TangoFilename   -- The raw output file generated by Tango 
       2) ProteinSequence -- The protein sequence analyzed by Tango. This sequence will have had all X's removed prior to analyses
       3) RemoveFile      -- By default, the program won't remove the original Tango file after executing. If the user wants
                             the Tango file to be deleted after analysis, use True when calling the function
**Output**: 

       All Tango scores with X's inserted in where the unknown residues were removed from the input protein in list format

### RunInterproscan.py

Author : Sara Willis
Date   : Wednesday February 27, 2019

This function will take in a fasta file containing protein sequences and will run them through InterProScan to determine the PfamIDs they are associated with as well as the expectation values

**INPUT**

     1) Fasta file containing the protein sequences to search against the database
     2) Path to InterProScan executable (interproscan.sh). This should be relative to the directory where this script is located
     3) Output Filename and path. This should be relative to the InterProScan directory, and should be directed to the 
     
**OUTPUT**

     Nested dictionary with the following format:
     {ProteinAccession: {PfamUID: [PfamUID1, ..., PfamUIDn], PfamStart: [PfamUID1_start, ..., PfamUIDn_start], PfamStop: [PfamUID1_stop, ..., PfamUIDn_stop], Eval: [PfamUID1_eval, ..., PfamUIDn_eval]}
