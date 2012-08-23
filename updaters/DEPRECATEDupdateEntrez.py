'''
Created on 17 Oct 2011

@author: Simon Bull
'''

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(entrezParsedOutput, schemaProteins, tableEntrezGene, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed GeneRIF data.
    #===========================================================================
    entrezData = utilities.file2list.main(entrezParsedOutput)
    entrezData = [eval(i) for i in entrezData]
    entrezDict = dict([(i[0], i) for i in entrezData])
    
    #===========================================================================
    # Extract the GeneRIF information recorded in the database.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableSELECT(cursor, '*', tableEntrezGene)
    results = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    
    #===========================================================================
    # Compare the parsed data with the data recorded in the table.
    #===========================================================================
    columnIndices = range(1, len(entrezData[0]))
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('SHOW COLUMNS FROM ' + tableEntrezGene)
    columns = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    columns = [i[0] for i in columns]
    
    toRemove = []
    toUpdate = {}
    toAdd = entrezDict.keys()
    for i in results:
        geneID = i[0]
        if entrezDict.has_key(geneID):
            # If the key is in both the parsed file and the table, then it does not need to be added.
            toAdd.remove(i[0])
            # Compare the row from the table with the parsed file, to determine if the table needs updating.
            for j in columnIndices:
                if i[j] != entrezDict[geneID][j]:
                    if not toUpdate.has_key(geneID):
                        toUpdate[geneID] = []
                    toUpdate[geneID].append(j)
        else:
            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
            toRemove.append(i[0])
    values = '(' + ('%s,' * len(entrezData[0]))
    values = values[:-1] + ')'
    
    #===========================================================================
    # Remove rows from the table that are not in the parsed file.
    #===========================================================================
    for i in toRemove:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.rowDELETE(cursor, tableEntrezGene, 'GeneID="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries removed from the Entrez Gene table: ', len(toRemove)
    
    #===========================================================================
    # Update rows that have different values in the parsed file and the table.
    #===========================================================================
    for i in toUpdate.keys():
        toSet = []
        for j in toUpdate[i]:
            updateString = columns[j] + ' = "' + entrezDict[i][j] + '"'
            toSet.append(updateString)
        toSet = ', '.join(toSet)
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableUPDATE(cursor, tableEntrezGene, toSet, 'GeneID="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries updated in the Entrez Gene table: ', len(toUpdate)
    
    #===========================================================================
    # Add rows which are not in the table, but are in the parsed file.
    #===========================================================================
    rowsToAdd = [entrezDict[i] for i in toAdd]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableINSERT(cursor, tableEntrezGene, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the Entrez Gene table: ', len(toAdd)