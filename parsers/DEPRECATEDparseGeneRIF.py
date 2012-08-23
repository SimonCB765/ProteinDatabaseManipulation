'''
Created on 10 Oct 2011

@author: Simon Bull
'''

def main(geneRIF):

    readIn = open(geneRIF, 'r')

    humanCount = 0
    humanGenes = set([])
    humanAnnotations = {}

    for line in readIn:
        if line[:4] == '9606':
            humanCount += 1
            chunks = line[:-1].split(None, 5)
            if humanAnnotations.has_key(chunks[1]):
                humanAnnotations[chunks[1]].append(chunks[-1])
            else:
                humanAnnotations[chunks[1]] = []
                humanAnnotations[chunks[1]].append(chunks[-1])
            humanGenes.add(chunks[1])

    readIn.close()

    return humanAnnotations