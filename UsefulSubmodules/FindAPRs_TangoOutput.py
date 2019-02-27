'''
Author : Sara Willis
Date   : Wednesday February 27, 2019

The purpose of this submodule is to search through a list of raw Tango scores to find runs of aggregation-prone amino acids. The function calculates the density of the AAs in APRs in the sequence to the length of the sequence, as well as the number of AAs in APRs in the sequence. These two values are returned to the user as well as a list of the specific runs. This is in case the user needs to calculate the aggregation scores for subsequences of the protein as well.

***** This differs from Scott's script:

If an unknown amino acid is found within an APR, the sequence is excluded from analyses and the user is returned False. If unknown amino acids exist anywhere else in the protein sequence, so long as they are not in a probable APR, the sequence is retained

Input: 
      1) AggregationList -- The raw Tango scores in list format. This can be obtained by using the module in the file ReplaceXInTangoOutput.py
      2) AggThreshold -- This is the value a Tango score must exceed to be counted as aggregation prone
      3) RunThreshold -- This is the number of Tango scores exceeding the AggThreshold in sequence for the region to be counted as an APR

Output:
     A list of lists. Each sublist contains the scores in every aggregation prone region.
'''

def FindAPRs(AggregationList, AggThreshold, RunThreshold):
    # runs is the master list that will hold each list belonging to an APR
    runs = []
    # We keep track of each individual potential APR with a temporary list
    run = []

    for aggValue in AggregationList:
        # a run is strung together from aggregation scores and X's
        if (type(aggValue) == float and float(aggValue) > AggThreshold) or aggValue == 'X':
            run.append(aggValue)
        # As soon as a score that doesn't exceed the AggTheshold is found, the run is broken
        else:
            # If the length of the run exceeds the RunThreshold, the run is checked to make sure there are no unknown AAs in the body of the run
            if len(run) >= RunThreshold:
                # Flanking unknown AAs are removed from the run and do not contribute to the length of the APR. As opposed to unknown AAs in the body of the run, these do not disqualify a sequence from analysis
                if run[0:2] == 'XX':
                    run = run[2:]
                if run[0] == 'X':
                    run = run[1:]
                if run[-2:] == 'XX':
                    run = run[:-2]
                if run[-1:] == 'X':
                    run = run[:-1]
                # After the flanking X's are removed from the run, if the length of the run exceeds the RunThreshold, then the sequence is checked to make sure there are no X's in the body of the sequence
                if len(run) >= RunThreshold:
                    if 'X' in run:
                        return False
                    else:
                        # If the run contains no unknown amino acids, then it is kept
                        runs.append(run)
                # If the run does not exceed the RunThreshold after the flanking unknown amino acids are removed, then the run is discarded and the aggregation list continues to be searched
                else:
                    pass
                runCount = 0
                run = []
    return runs


