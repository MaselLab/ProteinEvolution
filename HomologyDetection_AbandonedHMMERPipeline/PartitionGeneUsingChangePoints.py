'''

Author : Sara Willis
Date   : Wednesday February 27, 2019


This function is used to partition a gene based on change points in a multiple sequence alignment. It returns a dictionary of the segments of the original gene to the user.

The function takes in the ID of the sequence that needs to be partitioned, the change points and the alignment file.
'''

def PartitionOriginalGene(ChangePoints, MSA):
    # segmentNumber is used to differentiate the segments from one another
    # in the output dictionary
    newDictionary = {}

    # The multiple sequence alignment should be in fasta format 
    for record in SeqIO.parse('%s'%MSA, 'fasta'):
        segmentNumber = 0
        seq = str(record.seq) # sequence is gapped
        uid = record.id

        # Since the change points are relative to the alignment indices, the
        # function has to parse the gene relative to these coordinates

        # the list ChangePoints contains the indices of where the genes
        # should be partitioned. 
        for index in range(0,len(ChangePoints)-1):
            segmentNumber += 1
            # If a segment is not an empty string, it gets its own entry
            if seq[ChangePoints[index]:ChangePoints[index+1]].replace('-','') != '':
                newDictionary['%s_%s'%(uid,segmentNumber)] = seq[ChangePoints[index]:ChangePoints[index+1]].replace('-','')
    invalidKeys = []
    for key in newDictionary:
        # If it turns out a segment is less than 20 amino acids in length, then
        # an attempt is made to merge the segment with an adjacent segment. 
        if len(newDictionary[key]) <= 20:
            # a list of the keys that have sequence lengths less than 20 amino acids in length
            # are kept track of. This is used to remove UIDs later, since they cannot be removed
            # while a for loop is in progress
            invalidKeys.append(key)
            parsedKey = key.split('_')[0]
            segmentNumber = int(key.split('_')[1])
                
            try:
                if '%s_%s'%(parsedKey,segmentNumber+1) not in invalidKeys:
                    newDictionary['%s_%s'%(parsedKey,segmentNumber+1)] += newDictionary[key]
            except:
                # If the previous attempt to merge the sequence with another entry didn't
                # work, it's likely because the key segmentNumber +1 didn't exist (e.g., if
                # there are three segments and it's segment 3 that's too short, there's no
                # segment 4 to merge it with). If this is the case, then the segment is merged
                # with the preceding segment.
                try:
                    newDictionary['%s_%s'%(parsedKey,segmentNumber-1)] += newDictionary[key]
                # If both of the previous attempts to merge the small sequence with
                # another member of the dictionary have failed, then other keys don't
                # exist in the dictionary to merge the sequence with and it's left alone.

                # Note, this works even if there are two segments that have fewer than 20
                # amino acids since each entry is updated as it goes along, so nothing is lost
                except:
                    # If none of the above works, it's likely because the entire gene has
                    # fewer than 20 amino acids. If this is the case, the gene is kept as
                    # is
                    invalidKeys.remove(key)

    for key in invalidKeys:
        newDictionary.pop(key, None)
    dividedGenes = []
    for key in newDictionary:
        trimmedKey = key.split('_')[0].replace('R','').replace('L','')
        if trimmedKey not in dividedGenes:
            dividedGenes.append(trimmedKey)
    return newDictionary, dividedGenes
