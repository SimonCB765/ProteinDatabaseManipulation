'''
Created on 7 May 2011

@author: Simon Bull
'''

def main(toConvert, outputLoc):
    """Creates a FASTA format file from the information in toConvert.
    
    toConvert must be something like a list of lists, list of tuples, tuple of lists etc.
        toConvert[i][0] must be the UPID and toConvert[i][1] the sequence
    outputLoc is the location at which to write the FASTA file.
    
    """
    
    writeOut = open(outputLoc, 'w')
    
    for i in toConvert:
        writeOut.write('>' + str(i[0]) + '\n')
        writeOut.write(str(i[1]) + '\n')
    
    writeOut.close()
        