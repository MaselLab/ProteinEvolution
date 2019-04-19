import csv

'''
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


'''

def pull_tblout_hits(ECutoff, domECutoff, filename,bitBiasRatio=False):

    profileHits = set([])
    with open('%s'%filename, 'r') as f:
        reader = csv.reader(f, delimiter = ' ')
        for row in reader:
            # the space-delimiter used in the tblout can make for messy output, so we get rid of all the junk
            # and exclude any human-readable lines that don't contain data
            newRow = [i for i in row if i != '']
            if newRow[0][0] != '#':
                # Only a subset of the entries in tblout format are listed. More information on what's contained
                # in the output file can be found in the HMMER user manual.
                hit = newRow[0]
                fullEvalue = float(newRow[4])
                FullBit = float(newRow[5])
                FullBias = float(newRow[6])
                BestDomainEValue = float(newRow[7])
                bitBiasScore = FullBias/FullBit
                
                if BestDomainEValue > float(domECutoff) or fullEvalue >= float(ECutoff):
                    pass
                else:
                    if bitBiasRatio != False and bitBiasScore >= bitBiasRatio:
                        pass
                    else: 
                        profileHits.add(hit)
                    
    return profileHits
