

# The purpose of this function is to take in an input UID, the genes or gene segments
# that are hit when it is run through a homology detection pipeline, and a genes
# dictionary, and to run those through a change point detection program - See the fuction
# CreateBitScores and FindChangePoints above. It will use the output from the
# FindChangePoints function to determine the best change points to use.

# Input:
# inputUID - UID that will be partitioned
# results - list containing results from 3 iteration JackHMMER emulator run
# GenesDictionary - Relevant Dictionary to pull sequences from
def FindOptimalChangePoints(inputUID, results, GenesDictionary):
    fastaFilename = UserOptions()[7]
    msaFilename = UserOptions()[8]
    # a fasta file of the results is created
    for hit in results:
        CreateFastaFile(hit, GenesDictionary[hit], fastaFilename)
    # the hits are then aligned using clustal omega
    os.system('./clustalo --force --infile %s --outfile %s' %(fastaFilename, msaFilename))
    # and the fasta file is removed for decluttering
    os.remove(fastaFilename)

    # a bit score map is created using the alignment
    BitScores = CreateBitScores(msaFilename)
    # the change point penalty parameter and the minimum distance between change points
    # used with the FindChangePoints function are defined in the UserOptions module and
    # are loaded here.
    changePointPenalty = UserOptions()[16]
    minDistanceBetweenChangePoints = UserOptions()[18]

    # The change points are then found using the bit score map
    changePointOutput = FindChangePoints(BitScores, changePointPenalty, minDistanceBetweenChangePoints, UserOptions()[20])
    # The change points are a list and are the first entry in the function's output
    changePoints = changePointOutput[0]
    # and the average bit score values in each interval are the second entry in the output
    averages = changePointOutput[1]

    # If there are only two change points in the output, then that means there are no change points.
    # this is because the starting index (0) and the final index (len(bitscores)-1) are always included
    # as change points, so the minimum size of the change point list is 2. If only two entries
    # are included, then the gene doesn't need to be partitioned and the result False is returned
    # to the user.
    if len(changePoints) == 2:
        return False
    # if only one change point is found, then it is returned to the user. It is returned with True to
    # indicate that this was the only change point that was found.
    
    elif len(changePoints) == 3:
        return changePoints
    # if more than 2 change points are found, then we want to select two change points
    # based on the largest average in the bit score averages. This is because this is
    # usually associated with a region of homology.
    else:
        differences = [0]
        for changePoint in changePoints:
            if changePoint == 0:
                pass
            else:
                if changePoint == changePoints[-1]:
                    differences.append(0)
                    break
                else:
                    average1 = sum(BitScores[0:changePoint])/changePoint
                    average2 = sum(BitScores[changePoint:])/(BitScores[-1]-changePoint)
                    difference = abs(average2-average1)
                    differences.append(difference)
        maxDifference = max(differences)
        maxIndex = differences.index(maxDifference)
        bestChangePoint = changePoints[maxIndex]
        differences[maxIndex] = 0
        secondMaxDifference = max(differences)
        secondMaxIndex = differences.index(secondMaxDifference)
        secondBestChangePoint = changePoints[secondMaxIndex]
        parsedChangePoints1 = [0,bestChangePoint, changePoints[-1]]
        parsedChangePoints2 = [0, secondBestChangePoint, changePoints[-1]]
                              
        return parsedChangePoints1, parsedChangePoints2
