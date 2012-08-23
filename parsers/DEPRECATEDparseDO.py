'''
Created on 10 Oct 2011

@author: Simon Bull
'''

class Term(object):
    def __init__(self, termID, termName, synonyms, parentID):
        self.termID = termID
        self.termName = termName
        self.parentID = parentID
        self.synonyms = synonyms
        self.parent = None
        self.children = []

    def collect_children(self, childTerms=[]):
        childTerms.append(self.termName)
        for i in self.children:
            childTerms = i.collect_children(childTerms)
        return childTerms

    def collect_children_and_synonyms(self, childTerms=[]):
        childTerms.append(self.termName)
        childTerms.extend(self.synonyms)
        for i in self.children:
            childTerms = i.collect_children_and_synonyms(childTerms)
        return childTerms

    def tree_print_children(self, indent=''):
        print indent, self.termID, self.termName
        for i in self.children:
            i.tree_print_children(indent + '  ')

    def tree_print_parent(self):
        print self.termID, self.termName
        if self.parent != None:
            self.parent.tree_print_parent()

def parse(DOFile):

    termDict = {}  # Holds references to all the terms indexed by the term ID
    termDict['root'] = Term('root', 'root', 'root', 'root')

    ignore = True
    termID = ''
    termName = ''
    synonyms = []
    is_a = ''
    
    readFile = open(DOFile, 'r')

    lineCount = 0

    for line in readFile:
        lineCount += 1
        if line[:-1] == '[Term]':
            if not ignore:
                termDict[termID] = Term(termID, termName, synonyms, is_a)
            ignore = False
            termID = ''
            termName = ''
            synonyms = []
            is_a = 'root'
        elif line[:12] == 'is_obsolete:':
            ignore = True
        elif line[:8] == 'synonym:':
            chunks = line[:-1].split('"')
            synonyms.append(chunks[1])
        elif line[:3] == 'id:':
            chunks = line[:-1].split()
            termID = chunks[1]
        elif line[:5] == 'name:':
            chunks = line[:-1].split(None, 1)
            termName = chunks[1]
        elif line[:5] == 'is_a:':
            chunks = line.split()
            is_a = chunks[1]

    readFile.close()

    return termDict

def generate_tree(termDict):

    for termID in termDict:
        if termID != 'root':
            termDict[termID].parent = termDict[termDict[termID].parentID]
            termDict[termDict[termID].parentID].children.append(termDict[termID])

    return termDict

def main(DOFile):
    parseTree = parse(DOFile)
    parseTree = generate_tree(parseTree)
    
    return parseTree