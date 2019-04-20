# This function's purpose is to take in a dictionary of results from a HMMER search
# after the partitioning of the original gene, and to search to see if any of the
# results in the dictionary have a UID that occurs more than once, indicating a
# repeated element. If on is found, then True is returned to the user, otherwise
# False is returned if there are no repeated elements.
def RepeatedElementsIdentifier(resultsDictionary):
    newDictionary = {}
    try:
        for key in resultsDictionary:
            results = resultsDictionary[key]
            resultsWithoutAnnotation = [i.replace('L','').replace('R','') for i in results]
            resultsWithoutAnnotation = [i.split('_')[0] for i in resultsWithoutAnnotation]

            newDictionary[key] = resultsWithoutAnnotation
            for element in resultsWithoutAnnotation:
                if resultsWithoutAnnotation.count(element) > 1:
                    return True
        return False
    except:
        print('\n\n\nGenes not partitioned. Cannot identify repeated regions.\n\n\n')
        return False

