
'''
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
