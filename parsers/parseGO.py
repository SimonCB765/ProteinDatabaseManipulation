'''
Created on 7 Apr 2011

@author: Simon Bull
'''

import re

import utilities.MySQLaccess as mysql
import utilities.list2file

def main(goSchema, parsedGO):
    """
    Returns a file containing the parsed GO data.
    parsedGO - A file of 6-tuples, one on each line.
        The first element is the numerical identifier of the GO term.
        The second element is the name of the GO term.
        The third element is the category (biological_process, cellular component or molecular_function) that the term belongs to.
        The fourth element is all the paths from the term to the category it belongs to. The paths are separated from one another using semi-colons, and the elements of each path are separated from one another using '#'.
        The fifth element is all the level one terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the category in element three.
        The sixth element is all the level two terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the terms in element five.
    """

    # Connect to the GO schema.
    connGO, cursorGO = mysql.openConnection('root', goSchema)
    # Extract the information from the GO schema that is necessary to build up the paths.
    cursorGO = mysql.tableSELECT(cursorGO, 'id, name, term_type, acc', 'term')
    result = cursorGO.fetchall()

    # Turn the information extracted from the GO schema into a dictionary indexed by the identifiers of the GO terms.
    resultsDict = dict([(x[0], x[1:]) for x in result])

    # recordsToInsert will contain the information to insert into the GO path schema.
    recordsToInsert = []

    # For every term work out the name of the term, which of the three types of GO term it is, all the paths from the
    # term to the root terms (which are the three types), all unique level 1 and level 2 terms along the path.
    # A level 1 term is a child term of the root. A level 2 term is a child of a level 1 term.
    for r in resultsDict.keys():
        if re.match('GO:[0-9]+$', resultsDict[r][2]):
            # If the term_type is a GO term rather than something like a relationship describer (i.e. not something like
            # part_of, goslim_plant or regulates).
            pathDict = {}
            entry = resultsDict[r]
            startID = r
            startName = entry[0]
            startType = entry[1]
            pathDict[startID] = ['#' + startName]
            toCheck = [startID]

            while toCheck != []:
                # In order to determine all of the paths from a term to the root nodes, we treat the set of GO terms
                # as a directed graph. The direction of the edges is from parent to child. In the traversal we only go
                # in the direction opposite to that of the arrows (i.e. up a tree). Starting with a single term B, we find
                # the set of it's parent terms, P. To each element of P we add B to the end. The path for B then looks like
                # #P1#B #P2#B ... #Pn#B. This is repeated until all paths are capped at the front with a root term (the
                # root term will be the same for every path of B). Not all paths for B will be of the same length.
                current = toCheck.pop()
                currentPath = pathDict[current]
                # Get the terms that current is a subterm of.
                # Every entry in the term2term table that is returned represents one superterm of the current term.
                cursor = mysql.tableSELECT(cursorGO, 'term1_id', 'term2term', 'term2_id = \'' + str(current) + '\'')
                crossrefResults = cursor.fetchall()
                if crossrefResults == ():
                    # If the current term is not a subterm of any other term, then the current term is a root term.
                    # This means that the current path is finished.
                    rootNode = current
                    continue
                for i in crossrefResults:
                    # Get the id and name of the superterm being examined.
                    uniqueID = i[0]
                    name = resultsDict[uniqueID][0]
                    path = []
                    for j in currentPath:
                        # For every entry in the current path stick the superterm at the front of it, as it is an
                        # ancestor of whatever the head of the current path is.
                        path.append('#' + name + j)
                    # Get the identifiers of the terms which have already had their descendants enumerated.
                    keys = pathDict.keys()
                    if keys.count(uniqueID) == 0:
                        # If the current superterm has not already had its descendants enumerated, then simply record
                        # the paths just enumerated in the pathDict.
                        pathDict[uniqueID] = path
                    else:
                        # If the current superterm has already had its descendants enumerated, then add all the paths just
                        # enumerated to its record. Following this remove all duplicate paths.
                        pathDict[uniqueID].extend(path)
                        pathDict[uniqueID] = list(set(pathDict[uniqueID]))
                    # Add the superterm to the list of terms that need their superterms evaluating.
                    toCheck.append(uniqueID)

            # The final paths of interest are those recorded in the root node entry in the dictionary.
            final = pathDict[rootNode]
            levelOne = []
            levelTwo = []
            path = []

            for f in final:
                # For every path from the initial term to the root term, find the level 1 and level 2 terms.
                chunks = f.split('#')

                # Rather than setting the level 1 terms to be the 2nd and 3rd terms respectively, they are the 3rd and 4th.
                # This is because the GO database defines a term 'all' which is a dummy root node.
                if len(chunks) > 3:
                    levelOne.append(chunks[3])
                if len(chunks) > 4:
                    levelTwo.append(chunks[4])

                # Create the path that will be recorded in the database. Ignore the 'all' dummy root.
                path.append('#'.join(chunks[2:]))

            # Ensure only unique level 1 and level 2 terms are recorded.
            levelOne = list(set(levelOne))
            if levelOne == []:
                level1 = 'NA'
            else:
                level1 = ';'.join(levelOne)
            levelTwo = list(set(levelTwo))
            if levelTwo == []:
                level2 = 'NA'
            else:
                level2 = ';'.join(levelTwo)
            path = ';'.join(path)

            # Create the list of tuples to enter into the database.
            recordsToInsert.append((resultsDict[r][2][3:], startName, startType, path, level1, level2))

    mysql.closeConnection(connGO, cursorGO)

    utilities.list2file.main(recordsToInsert, parsedGO)