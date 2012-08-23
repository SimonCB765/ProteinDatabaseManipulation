'''
Created on 5 Dec 2011

@author: Simon Bull
'''

import utilities.MySQLaccess as mysql
import utilities.file2list

def main(pathwayElements, UPHumanNames, UPHumanAccessionMap, databasePassword, schemaProteins, tablePathways):
    
    # Generate the human accession mapping:
    UPAccMap = utilities.file2list.main(UPHumanAccessionMap)
    # Make a dictionary where the index is the, possibly deprecated UniProt accession, and the entry is the
    # representative accession.
    allUPID = {}
    for i in UPAccMap:
        chunks = i.split()
        allUPID[chunks[0]] = chunks[1]
    
    # Extract the names of the valid proteins.
    validProteinNames = utilities.file2list.main(UPHumanNames)

    # Extract the information about what pathways each element is in.
    proteinData = {}
    readIn = open(pathwayElements, 'r')
    for line in readIn:
        line = line.rstrip()
        chunks = line.split('\t')
        pathway = chunks[0]
        elements = chunks[2:]
        for i in elements:
            # For every element in the pathway, extract the names and accessions of the protein. Also record the
            # name of the pathways that the protein is associated with.
            chunks = i.split(':')
            cpath = chunks[0]
            type = chunks[1]
            name = chunks[2]
            UPAccession = chunks[3]
            if proteinData.has_key(cpath):
                proteinData[cpath]['Names'] |= set(name.split(','))
                proteinData[cpath]['UPAccessions'] |= set(UPAccession.split(','))
                proteinData[cpath]['Pathways'].add(pathway)
            else:
                proteinData[cpath] = {'Names' : set(name.split(',')), 'UPAccessions' : set(UPAccession.split(',')), 'Pathways' : set([pathway])}
    readIn.close()
    
    # Consolidate the pathway element information extracted. This means combining records with UniProt accessions
    # that map to the same accession, and removing any proteins that can't be found to have a known accession or
    # a valid name.
    consolidatedProteinData = {}
    for i in proteinData.keys():
#        names = proteinData[i]['Names']
        UPAccessions = proteinData[i]['UPAccessions']
        pathways = proteinData[i]['Pathways']
        for j in UPAccessions:
            j = j.split('-')[0]  # Combine all isoforms into one record.
            # Convert each UniProt accession into its representative form.
            if allUPID.has_key(j):
                reprAcc = allUPID[j]
                if consolidatedProteinData.has_key(reprAcc):
                    consolidatedProteinData[reprAcc]['UPAccessions'] |= UPAccessions
                    consolidatedProteinData[reprAcc]['Pathways'] |= pathways
                else:
                    consolidatedProteinData[reprAcc] = {'UPAccessions' : UPAccessions, 'Pathways' : pathways}
            else:
                # Don't use the names for anything right now.
                pass
    
    proteinTuples = []
    for i in consolidatedProteinData.keys():
        proteinTuples.append(tuple([i, len(consolidatedProteinData[i]['Pathways'])]))
    
    values = '(' + ('%s,' * len(proteinTuples[0]))
    values = values[:-1] + ')'
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    cursor.execute('TRUNCATE TABLE ' + tablePathways)
    mysql.tableINSERT(cursor, tablePathways, values, proteinTuples)
    mysql.closeConnection(conn, cursor)