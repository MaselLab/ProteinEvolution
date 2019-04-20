import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


'''
The purpose of this function is to create a fasta file with sequences designated by the user. The function taskes in three arguments:

 1) The identifier for the sequence in the fasta file
 2) Protein sequence
 3) The saved filename

 IMPORTANT:

The function writes to the file using the "append" function, so sequence can be written to the file by calling the function in a for loop. This, however, means that the user must remember to delete old files as, if they are not, new data may be inadvertently appended to the ends of old files 
'''

def CreateFastaFile(UID, Protein, HitsFastaFile):
    # we save sys.stdout as a new variable so that we can redirect what's printed
    # to the terminal to a file, and then revert to printing to the terminal at the
    # end of the function's use
    orig_stdout = sys.stdout
    record = SeqRecord(Seq(Protein), id= '%s'%UID)
    # stdout is redirected to a file
    sys.stdout = open('%s'%HitsFastaFile, 'a')
    print(record.format('fasta'))
    # Once it has been printed, sys.stdout is closed and the program will start
    # printing to the terminal again
    sys.stdout.close()
    sys.stdout=orig_stdout
