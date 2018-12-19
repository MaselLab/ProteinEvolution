
# ExtractNCBIGenomes_UploadResultsToMySQL

----------------------------------------------------------------------

FullPipeline.py is a command-line based python script that does the following:

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

----------------------------------------------------------------------

## Word of Caution!

Once this script starts running in full with the supplied species list, it takes approximately three weeks to finish. Be *really* sure that everything in the pipeline is functioning before starting. 

There is an option on line 156 that allows the user to run only ten genes per species to ensure each step of the pipeline is functioning before a full run. This should be removed before executing the script for entire genomes.

----------------------------------------------------------------------

## Dependencies: 

	Local installation of InterProScan

	MySQL

	IUPred* (included)
		-- The files hist, histo_casp, hist_sum, ss, and ss_casp are required for iupred to function

*The Pairwise Energy Content Estimated from Amino Acid Composition Discriminates between Folded and Intrinsically Unstructured Proteins 
Zsuzsanna Dosztányi, Veronika Csizmók, Péter Tompa and István Simon 
J. Mol. Biol. (2005) 347, 827-839.


----------------------------------------------------------------------
Python modules required:

	mysql.connector, Biopython
