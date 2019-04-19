import math
from Bio import SeqIO

'''
The purpose of this submodule is to find the heights of the sequence logo of a multiple sequence alignment.

The function accepts a multiple sequence alignment in fasta format as input. Each column in the alignment is used to generate the sequence logo array. This submodule can be used with a change point detection module to find optimal regions to partition a multiple sequence alignment.

For more information on sequence logos and the mathematics used to generate the map: https://en.wikipedia.org/wiki/Sequence_logo

Note: The bit score, when plotted, should look identical to the heights of the sequence logo map that shows up in Geneious for multiple sequence alignments
'''
def CreateBitScores(InputMSAFile):
    # The start variable tells the function that it is looking at the first residue
    # in a column
    start = True
    # The columns of the MSA are stored as sublists in a master list. The master
    # list is declared here
    Columns = []
    # N keeps track of the number of sequences in the alignment
    N = 0
    # The MSA is parsed using the BioPython module
    for record in SeqIO.parse('%s'%InputMSAFile, 'fasta'):
        N += 1
        UID = record.id
        sequence = record.seq
        # If the sequence is the first to be added to the column, a list for that
        # column is created and added to the overall list
        if start == True:
            # Each residue in the sequence is checked
            for n in range(0,len(sequence)):
                # if the sequence is gapped for a specific column, nothing is added
                # to the column list and an empty list is added to the column position
                # in the overall list
                if sequence[n] == '-':
                    Columns.append([])
                else:
                    # if it's ungapped for the specific column, then a list with the
                    # relevant residue is added to the overall list.
                    Columns.append([sequence[n]])
            # once the initial sequence in the alignment is parsed, start is set to false and
            # all subsequent residues will be appended to the preexisting column sublists.
            start = False
        else:
            for n in range(0,len(sequence)):
                if sequence[n] == '-':
                    pass
                else:
                    Columns[n].append(sequence[n])                
    # The error correction score is the following
    e_n = (1/math.log(2))*(19/(2*N))
    # each column will have a bit score/height associated with it that represents the height of the
    # overall sequence logo. If you would like more information on the mathematics, visit the Wikipedia
    # page referenced in this script's description
    BitScores = []
    # Each column in the MSA has a bit score assigned to it, that must be calculated based on the residues
    # that occur in that column
    for column in Columns:
        # H_i is the uncertainty, or Shannon entropy, of position i
        H_i = 0
        bitScore = 0
        # We only want to search for residues that exist for the column, so we assemble a unique list
        Residues = []
        for residue in column:
            if residue not in Residues:
                Residues.append(residue)
        for residue in Residues:
            # For each residue that shows up, we want to know the number of times it occurs in that column
            residueCount = column.count(residue)
            # We get the normalized count by dividing the number of occurences by the number of sequences in
            # the MSA
            normalizedCount = residueCount/N
            H_i -= (normalizedCount * math.log(normalizedCount,2))
        # We then calculate the information content of the column using all uncertainty scores and the error
        # error correction score. Log_2(20) is specifically for amino acid alignments. This would be
        # Log_2(4) for nucleotides
        R_i =  math.log(20,2) - (H_i + e_n)
        for residue in Residues:
            residueCount = column.count(residue)
            normalizedCount = residueCount/N
            # The bit score is then calculated as the sum of all normalized counts (for each occuring residue
            # in a column) times the information content we calculated.
            bitScore += normalizedCount*R_i
        # The bit scores are then stored in a list with one entry per column in the MSA
        BitScores.append(bitScore)
    return BitScores
