'''
Created on 18 Apr 2011

@author: Simon Bull
'''


def main(fileLocation):
    """Parses a default Pepstats output file, and returns information as a dictionary.
    
    The output file must be generated from a FASTA format file that has any number of sequences in it.
    
    @param fileLocation: the location of the Pepstats output file
    @type fileLocation:  string
    
    """
    
    resultDict = {}
    
    readResults = open(fileLocation, 'r')
    
    for line in readResults:
        
        if line[:11] == 'PEPSTATS of':
            # Got the start of a new protein results. If a protein is in the results file multiple times
            # then resultDict will only contain one instance of it.
            chunks = line.split()
            currentProt = chunks[2]
            resultDict[currentProt] = {}
        elif line[:17] == 'Isoelectric Point':
            # If this is true the isoelectric point information has been found.
            chunks = line.split()
            pI = chunks[3]
            resultDict[currentProt]['pI'] = pI
    
    readResults.close()
    
    return resultDict