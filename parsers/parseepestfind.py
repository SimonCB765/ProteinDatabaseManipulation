'''
Created on 18 Apr 2011

@author: Simon Bull
'''

def main(fileLocation):
    """Parses a default epestfind output file, and returns the number of valid PEST motifs found.
    
    @param fileLocation: the location of the epestfind output file
    @type fileLocation:  string
    
    """
    
    numberOfPEST = 0
    
    readResults = open(fileLocation, 'r')
    
    for line in readResults:
        if line[:20] == 'Potential PEST motif':
            numberOfPEST += 1
    
    readResults.close()
    
    return numberOfPEST