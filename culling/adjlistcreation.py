'''
Created on 5 Feb 2011

@author: Simon Bull
'''

import sparsematrix

def main(similarities, cutoffPercent=20, maxEValue=1, minAlignLength=20):
    """Create a sparse matrix from the processed PSI-BLAST output.
    
    @param similarities: The location of the processed PSI-BLAST output.
    @type similarities: string
    @param cutoffPercent: A percentage similarity > this parameter is deemed to be too similar.
    @type cutoffPercent:  float
    @param maxEValue: The maximum permissible value for the E value of an alignment.
    @type maxEValue: float
    @param minAlignLength: The number of amino acids aligned in the query and the hit sequence
                           must be >= this value for the percentage similarity to be deemed significant.
    @type minAlignLength: integer
    
    """

    proteinNames = []  # Store the names of all the proteins found to be too similar to another protein
    similarProteins = []  # Store the pairs that are too similar
    
    readSimilarities = open(similarities, 'r')
    
    for line in readSimilarities:
        chunks = line.split()
        query = chunks[0]
        hit = chunks[1]
        percentage = chunks[2]
        alignmentLength = chunks[3]
        EValue = chunks[4]
        
        # Ignore similarities where the query and the hit are the same, the percentage similarity is <= cutoffPercent,
        # the E value is > maxEValue and the length of the alignment is < minAlignLength.
        invalid = (query == hit or
                   float(percentage) <= cutoffPercent or
                   float(EValue) > float(maxEValue) or
                   int(alignmentLength) < minAlignLength
                   )

        if not invalid:
            # If the similarity is valid record the proteins as being too similar.
            proteinNames.append(query)
            proteinNames.append(hit)
            similarProteins.append(tuple(sorted((query, hit))))
    
    readSimilarities.close()
    
    proteinNames = list(set(proteinNames))
    proteinNames.sort()
    similarProteins = list(set(similarProteins))
    indexDict = dict((proteinNames[x], x) for x in range(len(proteinNames)))
    
    # Create the sparse matrix
    adjacent = sparsematrix.SparseMatrix(len(proteinNames))
    xValues = [indexDict[x] for (x,y) in similarProteins]
    yValues = [indexDict[y] for (x,y) in similarProteins]
    adjacent.addlist(xValues, yValues)
    adjacent.addlist(yValues, xValues)
    
    return adjacent, proteinNames
