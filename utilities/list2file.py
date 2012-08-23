'''
Created on 6 Apr 2011

@author: Simon Bull
'''

def main(toWrite, outputFile):
    """Outputs the list toWrite to the file outputFile with one element on each line.

    toWrite @type - list
    toWrite @use  - the list to be output to file with one element of the list on each line
    outputFile @type - string
    outputFile @use  - the pathname of the location to write the list out to
    """

    if type(toWrite) != list:
        raise TypeError("The first argument must be a list.")

    writeTo = open(outputFile,'w')

    for i in toWrite:
        print>>writeTo, i
    writeTo.close()