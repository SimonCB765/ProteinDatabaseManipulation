'''
Created on 6 Apr 2011

@author: Simon Bull
'''

def main(inputFile):
    """Used to convert a file of elements into a list.

    Note: the file should have one element of the list on each line.

    inputFile @type - string
    inputFile @use  - location of the file to convert to a list
    return @type - list
    return @use  - each element will be one line of the file with the trailing '\n' stripped off
    """

    readFrom = open(inputFile,'r')
    outputList = []

    for i in readFrom:
        if i[-1] == '\n': #Strip out \n
            outputList.append(i[:-1])
        else:
            outputList.append(i)

    readFrom.close()
    return outputList