import ruptures as rpt
import matplotlib.pyplot as plt


'''
The purpose of this function is to read in a list of numerical values and output the detected change points associated with that list. For example, if the user wanted to know how to partition a multiple sequence alignment, they could create a bit score list, feed it into this module, and determine where the alignment should be split. The function has three inputs that need to be specified by the user, as well as one optional parameter:

 [1] The bit scores list is the first argument in the function and is what is analyzed for change points
 [2] penalty is what is used to determine how sensitive the program should be when looking for change points. Too sensitive means every slight variation is counted as a change point, while very strict means some may be missed and only the most extreme changes are counted, so some fiddling may be necessary
 [3] minimumSize is the minimum number of elements allowed between change points. This can prevent hundreds of change points showing up in a small bumpy region
 [4] The last argument is optional. The user can opt to set this to True if they wish to see a display of the bit scores along with the change points. This can be helpful during testing/debugging so the user can better visualize the data points. 

 The output of this function is a list of the coordinates where change points were found. The start and stop locations in the MSA are always counted as change points for convinience. A second list of the averages of the bit scores in each partitioned section of the MSA, divided up by the change points, is also returned to the user. This makes spikes easier to identify. It also allows the user to pick the best change point (most dramatic difference between regions), if the user only wants to focus on one region at a time.

Dependencies:

Ruptures: a change point detection module in Python. 
          -- http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html
          -- https://pypi.org/project/ruptures/

'''

def FindChangePoints(BitScores, penalty, minimumSize, Display=False):
    # This function uses the module Ruptures to detect change points in the bit score map.
    algo = rpt.Pelt(min_size = minimumSize).fit(np.array(BitScores))
    result = algo.predict(pen = penalty)
    # After the change points are found using ruptures, the program makes sure that the last
    # index in the MSA is included as a change point
    if len(BitScores)-1 not in result:
        result.append(len(BitScores)-1)
    # Because python starts its counting at 0, if the length of the MSA is included as
    # a result, it will cause an indexing error if it's used. This is why it's removed
    # if it exists.
    if len(BitScores) in result:
        result.remove(len(BitScores))
    # It also makes sure the starting index is in the results
    if 0 not in result:
        result.append(0)
    result = sorted(result)

    # The average bit score value in each partitioned section of the bit score map
    # is then determined
    averages = []
    for changePoint in result:
        # If the last change point is found, then there will not be a change point after it,
        # so the loop ends
        if result.index(changePoint) == len(result)-1:
            break
        # otherwise, the change point after the change point selected is chosen 
        else:
            changePoint2 = result[result.index(changePoint)+1]
        # The total bit score is set to zero
        total = 0
        # and then each bit score between the two change points is added to the total
        for score in BitScores[changePoint:changePoint2]:
            total += score
        # The average is then found by dividing the total by the length of the interval
        average = total/(changePoint2-changePoint)
        # and the average is added to the list that will be returned to the user
        averages.append(average)
    # If the user has specified that they would like to see the change points
    # illustrated on the change point map, then the plot is printed.
    if Display == True:
        rpt.display(np.array(BitScores), result)
        plt.show()

    # The program then returns the list of change points and the average bit score values
    # between them to the user.
    return result, averages
