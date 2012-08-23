'''
Created on 19 Apr 2011

@author: Simon Bull
'''

def main(outputSEG):
    """Used to parse the output of the SEG masking program to determine the number of low complexity regions.
    
    The output must be generated from a FASTA format file of a single protein sequence. The output format of the SEG
    masking program must be the interval format. The interval format output generated froma FASTA file is as follows:
    >sequence
    startInterval1 - endInterval1
    startInterval2 - endInterval2
    ...
    startIntervalN - endIntervalN
    
    @param outputSEG: the location of the SEG output file
    @type outputSEG:  string
    
    """
    
    numRegions = 0
    
    readIn = open(outputSEG, 'r')
    
    for line in readIn:
        if line[0] == '>':
            continue
        elif len(line) > 1:
            # If the line is not empty
            numRegions += 1
    
    return numRegions