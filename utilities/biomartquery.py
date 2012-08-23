'''
Created on 13 Oct 2011

@author: Simon Bull
'''

import subprocess
import os

def single_filter(martName, datasetName, filterName, filterValue, attributeNames, biomartQueryScript, bioMartXMLQuery,
                  outputFile):
    
    # Delete the file if it exists.
    if os.path.isfile(outputFile):
        os.remove(outputFile)
    

    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query client=\"pythonclient\" processor=\"TSV\" limit=\"-1\" header=\"0\" uniqueRows = \"1\">
\t<Dataset name=\"%s\" config=\"%s\">
\t\t<Filter name=\"%s\" value=\"%s\" />
""" % (datasetName, martName, filterName, filterValue)
    
    for i in attributeNames:
        query += '\t\t<Attribute name=\"%s\" />\n' % (i)

    query = query[:-1] + """
\t</Dataset>
</Query>"""

    writeOut = open(bioMartXMLQuery, 'w')
    for line in query:
        writeOut.write(line)
    writeOut.close()

    subprocess.call('perl ' + biomartQueryScript + ' ' + bioMartXMLQuery + ' ' + outputFile)

def multi_filter(martName, datasetName, filters, attributeNames, biomartQueryScript, bioMartXMLQuery, outputFile):
    
    # Delete the file if it exists.
    if os.path.isfile(outputFile):
        os.remove(outputFile)
    

    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query client=\"pythonclient\" processor=\"TSV\" limit=\"-1\" header=\"0\" uniqueRows = \"1\">
\t<Dataset name=\"%s\" config=\"%s\">
""" % (datasetName, martName)

    for i in filters.keys():
        query += '\t\t<Filter name=\"%s\" value=\"%s\" />\n' % (i, filters[i])
    
    for i in attributeNames:
        query += '\t\t<Attribute name=\"%s\" />\n' % (i)

    query = query[:-1] + """
\t</Dataset>
</Query>"""

    writeOut = open(bioMartXMLQuery, 'w')
    for line in query:
        writeOut.write(line)
    writeOut.close()

    subprocess.call('perl ' + biomartQueryScript + ' ' + bioMartXMLQuery + ' ' + outputFile)