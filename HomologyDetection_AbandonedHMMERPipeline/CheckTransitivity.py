# The purpose of this function is to determine if one set of results is a subset of the other. This will determine whether an alignment is split or not
def Transitivity(segmentedResults):
    Transitivity = True
    # If there are only two sets of results in the SegmentedResults dictionary,
    # then the input alignment was only partitioned into two sections. If this
    # is the case, then comparing the results is easy
    if len(segmentedResults) == 2:
        first = True
        for result in segmentedResults:
            if first == True:
                result1 = [i.split('_')[0].replace('R','').replace('L','') for i in segmentedResults['1']]
                first = False
            else:
                result2 = [i.split('_')[0] for i in segmentedResults[result]]
        if set(result1).issubset(set(result2)) == True or set(result2).issubset(set(result1)) == True:
            Transitivity = True
            return Transitivity
        else:
            Transitivity = False
            return Transitivity
    # If there are more than two sets of results in the SegmentedResults dictionary,
    # then there are three (there can't be more given how the algorithm works). This
    # means there are three pairwise comparisons that need to be made
    else:
        result1 = [i.split('_')[0].replace('R','').replace('L','') for i in segmentedResults['1']]
        result2 = [i.split('_')[0].replace('R','').replace('L','') for i in segmentedResults['2']]
        result3 = [i.split('_')[0].replace('R','').replace('L','') for i in segmentedResults['3']]

        if set(result1).issubset(set(result3)) == False and set(result3).issubset(set(result1)) == False:
            Transitivity = False
            return Transitivity
        elif set(result1).issubset(set(result2)) == False and set(result2).issubset(set(result1)) == False:
            Transitivity = False
            return Transitivity
        elif set(result2).issubset(set(result3)) == False and set(result3).issubset(set(result2)) == False:
            Trasitivity = False
            return Transitivity
        else:
            Transitivity = True
            return Transitivity
                                   
