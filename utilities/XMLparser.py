'''
Created on 16 Sep 2011

@author: Simon Bull
'''

import re

class Tag(object):
    """A class to record the tags in an xml file.

    Each Tag object functions as a node in a tree. It contains information about the parent Tag, children Tags and
    information about the parameters found in the tag heading and the data potentially recorded in the tag.
    """

    def __init__(self, tagName, parent, parameters=''):
        self.tagName = tagName
        self.parent = parent
        self.children = []
        self.parameters = {}
        if parameters != '':
            chunks = re.findall('\S+=".*?"', parameters)  # Matches a positive number of non-whitespace, then =, then ", then the smallest amount of
                                                          # characters before the next ", then the "
            for i in chunks:
                i = i.replace('"', '')
                splits = i.split('=')
                self.parameters[splits[0]] = splits[1]
        self.data = ''

    def add_child(self, child):
        """Add a child (descendant) to the list of children."""
        self.children.append(child)

    def add_data(self, data):
        """Append data to the currently held string of data."""
        self.data += str(data)

    def tree_print_all(self, indent=''):
        """Print the tags, parameters and data for the tree starting at the current Tag."""
        print indent, self.tagName,
        for i in self.parameters:
            print i, '=', self.parameters[i], ' ',
        print ''
        if self.data != '':
            print indent, '  ', repr(self.data)
        for i in self.children:
            i.tree_print_all(indent + '  ')

    def tree_print_tags(self, indent=''):
        """Print the tags for the tree starting at the current Tag."""
        print indent, self.tagName
        for i in self.children:
            i.tree_print_tags(indent + '  ')

    def tree_print_tags_data(self, indent=''):
        """Print the tags and data for the tree starting at the current Tag."""
        print indent, self.tagName
        if self.data != '':
            print indent, '  ', repr(self.data)
        for i in self.children:
            i.tree_print_tags_data(indent + '  ')

    def tree_print_tags_params(self, indent=''):
        """Print the tags and parameters for the tree starting at the current Tag."""
        print indent, self.tagName,
        for i in self.parameters:
            print i, '=', self.parameters[i], ' ',
        print ''
        for i in self.children:
            i.tree_print_tags_params(indent + '  ')


class XMLParser(object):
    """A class to contain a xml file, along with the results of its parsing."""
    def __init__(self, XMLFile):
        self.XMLFile = XMLFile
        self.parseTree = None
        self.accessTree = None

    def load_file(self, XMLFile):
        """Load a file into the XMLParser object."""
        self.XMLFile = XMLFile
        self.parseTree = None
        self.accessTree = None

    def match_replacement(self, matchobj):
        """Used to replace matches when parsing the XML tags."""
        return 'a' * len(matchobj.group(0))

    def parse(self):
        """Parse the xml file recorded in the XMLParser object."""
        self.parseTree = Tag('dummy_root', None)
        self.accessTree = {}

        readXML = open(self.XMLFile, 'r')

        tagStack = [self.parseTree]
        accessStackName = []
        accessStackChildPos = []
        currentParent = None
        currentTag = self.parseTree
        data = ''

        # The parsing works by going through each line of the xml file character by character. The names of tags in the current
        # path are recorded in a stack of tag names (i.e. if the file is:
        # <drugs>
        #     <drug>
        #         <name>
        # the current stack would be ['drugs', 'drug'] if we are about to look at <name>.
        # When a new tag is found a new Tag object is created. This new tag object is indexed within the parse tree by
        # concatenating the elements of the stack together with a '.' in between each (e.g. here the index for the name tag
        # would be 'drugs.drug.name'). This entry ('drugs.drug.name') is placed in an access dictionary which enables easy
        # access to the Tag objects in the parse tree. The access tree contains for each Tag object in the parse tree, the
        # path that can be followed from the root to reach it. In this example 'drugs.drug.name' can be reached by:
        # parseTree.children[0].children[0]
        # Therefore, the entry in the eaccess tree for 'drugs.drug.name' is [0,0].

        for line in readXML:
            if line[:2] == '<?':
                # Skip the definition tags that start with a <?
                continue
            if line[-1] == '\n':
                line = line[:-1]  # Strip the newline character.
            line = line.lstrip()  # Remove the left tabs and spaces used in XML.
            currentPos = 0
            maxPosition = len(line)
            if data != '':
                # If there is data for the current tag then record it.
                data += ' '
            # Read through the line in the xml file recording the start and end of tag blocks, the parameters and the data
            # for the tags.
            while currentPos < maxPosition:
                if line[currentPos] == '<':
                    # If the line pertains to a tag, either beginning or ending one.
                    if line[currentPos + 1] == '/':
                        # Found the end of a tag.
                        tagEnd = re.match('</.*?>', line[currentPos:])
                        tagEnd = tagEnd.group(0)
                        # Add (potentially) any data associated with the tag.
                        data = data.strip()
                        currentTag.add_data(data)
                        data = ''
                        tagStack.pop()  # Pop the current tag off the tag stack as it's no longer needed.
                        currentParent.add_child(currentTag)
                        currentTag = currentParent
                        if len(tagStack) == 1:
                            # Only happens when closing the final tag, so that the dummy root is the only tag left on the stack.
                            currentParent = None
                        else:
                            currentParent = tagStack[-2]

                        accessStackName.pop()
                        accessStackChildPos.pop()

                        currentPos += len(tagEnd)
                    else:
                        # Found the start of a tag, or an empty tag.
                        # When determining the start and end position of a tag, it's necessary to ignore all
                        # characters in parameter definitions (i.e. inside the ""). This is because the parameter
                        # definitions may contain reserved characters, specifically they may contain a >, and falsely
                        # indicate the end of a tag.
                        tagStart = re.match('<.*?>', re.sub('".*?"', self.match_replacement, line[currentPos:]))
                        tagStart = line[currentPos : currentPos + len(tagStart.group(0))]
#                        tagStart = tagStart.group(0)
                        chunks = tagStart[1:-1].split(None, 1)
                        tagName = chunks[0]
                        currentParent = currentTag
                        # Determine if the tag is empty.
                        if tagStart[-2:] == '/>' and chunks == 0:
                            # Found an empty element tag, which should be ignored. The tag is empty if it opens and
                            # closes on one line, and it has no parameters e.g. <A/>
                            currentPos += len(tagStart)
                            continue
                        elif tagStart[-2:] == '/>':
                            # This is true if the tag opens and closes on the same line, and has parameters.
                            # Work around this by saving the params, and then pretending the end of a tag has
                            # been found. For example, if the line is <A id="10"/>, then record the tag as normal.
                            # However, need to alter the line to be <A id="10"</>, and also need to set the parser's
                            # pointer to the beginning of the </>. This will allow the nest iteration of the while loop
                            # to start on the </>, and recgnise that a tag end has been reached.
                            tagStart = 'a' * (len(tagStart) - 2)
                            line = line[:-2] + '</>'
                        if len(chunks) > 1:
                            # If this is true then the tag contains parameters in its definition.
                            currentTag = Tag(tagName, currentParent, parameters=chunks[1])
                        else:
                            currentTag = Tag(tagName, currentParent, parameters='')
                        tagStack.append(currentTag)

                        accessStackName.append(tagName)
                        accessStackChildPos.append(len(currentParent.children))  # Get the location that the newly found Tag
                                                                                 # will take in the record of children of the
                                                                                 # parent Tag.
                        if not self.accessTree.has_key('.'.join(accessStackName)):
                            # If the tree recording the structure of the parse tree does not have an entry for the current
                            # tag then add it.
                            self.accessTree['.'.join(accessStackName)] = [list(accessStackChildPos[1:])]
                        else:
                            self.accessTree['.'.join(accessStackName)].append(list(accessStackChildPos[1:]))

                        currentPos += len(tagStart)
                else:
                    # Not reading a tag, so must be reading data.
                    data += line[currentPos]
                    currentPos += 1

        readXML.close()

        self.parseTree = self.parseTree.children[0]

    def retrieve(self, terms):
        """Retrieves the Tag objects associated with the terms desired.

        Terms take the form X.Y.Z where X, Y and Z are tags in the XML file, and Z is a child of Y which is a child of X."""

        # Extract the access paths for the terms of interest.
        accessPaths = dict([(i, self.accessTree[i]) for i in terms])

        # Determine the number of different data entries
        dataEntries = range(len(set([i[0] for i in accessPaths[terms[0]]])))
        dataToReturn = dict([(i, dict([(j, []) for j in terms])) for i in dataEntries])

        for i in accessPaths:
            traversal = accessPaths[i]
            for j in traversal:
                # Determine which data entry the traversal path corresponds to based on the first number in the traversal list.
                dataEntry = j[0]
                current = self.parseTree
                for k in j:
                    current = current.children[k]
                dataToReturn[dataEntry][i].append(current)

        return dataToReturn

    def retrieve_data(self, terms):
        """Retrieves the data associated with the tags corresponding to the terms desired.

        Terms take the form X.Y.Z where X, Y and Z are tags in the XML file, and Z is a child of Y which is a child of X.
        If one of the terms has no data, then it is returned as ......."""

        tags = self.retrieve(terms)

        dataToReturn = {}

        for i in tags:
            dataToReturn[i] = {}
            for j in tags[i]:
                dataToReturn[i][j] = {}
                dataToReturn[i][j]['data'] = []
                for k in tags[i][j]:
                    dataToReturn[i][j]['data'].append(k.data)

        return dataToReturn

    def retrieve_data_and_params(self, terms):
        """Retrieves the data and parameters associated with the tags corresponding to the terms desired.

        Terms take the form X.Y.Z where X, Y and Z are tags in the XML file, and Z is a child of Y which is a child of X.
        If one of the terms has no data, then it is returned as .......
        If one of the terms has no parameters the value for its parameters is returned as an empty dictionary."""

        tags = self.retrieve(terms)

        dataToReturn = {}

        for i in tags:
            dataToReturn[i] = {}
            for j in tags[i]:
                dataToReturn[i][j] = {}
                dataToReturn[i][j]['data'] = []
                dataToReturn[i][j]['params'] = []
                for k in tags[i][j]:
                    dataToReturn[i][j]['data'].append(k.data)
                    dataToReturn[i][j]['params'].append(k.parameters)

        return dataToReturn

    def retrieve_params(self, terms):
        """Retrieves the parameters associated with the tags corresponding to the terms desired.

        Terms take the form X.Y.Z where X, Y and Z are tags in the XML file, and Z is a child of Y which is a child of X.
        If one of the terms has no parameters the value for its parameters is returned as an empty dictionary."""

        tags = self.retrieve(terms)

        dataToReturn = {}

        for i in tags:
            dataToReturn[i] = {}
            for j in tags[i]:
                dataToReturn[i][j] = {}
                dataToReturn[i][j]['params'] = []
                for k in tags[i][j]:
                    dataToReturn[i][j]['params'].append(k.parameters)

        return dataToReturn