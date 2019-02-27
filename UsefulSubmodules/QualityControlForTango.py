'''

Author : Sara Willis
Date   : February 27, 2019


This function is used to check to see whether a sequence passes imposed quality filters before it's run through the aggregation-prediction software Tango. Specifically, if a protein has three or more unknown residues in sequence, then the protein does not pass the filter and the user is notified to not use the sequence for Tango analysis. 

Input: A protein sequence

Output: True or False depending on whether the sequence passes the quality check or not.
'''

def PassesTangoQualityControl(Sequence):

    # It will look for runs of X's (unknown residues) of three or more in length.
    XRun = 0
    for residue in Sequence:
        if residue == 'X':
            XRun += 1
        else:
            if XRun >= 3:
                return False
            else:
                XRun = 0
    # If no runs >= 3 were found, the sequence passes the quality check and the user is returned True             
    return True
