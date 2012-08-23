'''
Created on 18 Oct 2011

@author: Simon Bull
'''

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(unigeneParsedOutput, unigeneParsedTotals, schemaProteins, tableUniGene, tableUniGeneTotals, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed UniGene data, and the UniGene expression totals.
    #===========================================================================
    UGData = utilities.file2list.main(unigeneParsedOutput)
    UGData = [tuple([int(j) for j in eval(i)]) for i in UGData]
    UGDict = dict([(i[0], i) for i in UGData])
    
    UGTotalsData = utilities.file2list.main(unigeneParsedTotals)
    UGTotalsData = [tuple([j[0], eval(j[1])]) for j in [eval(i) for i in UGTotalsData]]
    
    #===========================================================================
    # Extract the UniGene information recorded in the database.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableSELECT(cursor, '*', tableUniGene)
    results = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    
    #===========================================================================
    # Compare the parsed data with the data recorded in the expression table.
    #===========================================================================
    columnIndices = range(1, len(UGData[0]))
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('SHOW COLUMNS FROM ' + tableUniGene)
    columns = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    columns = [i[0] for i in columns]
    
    toRemove = []
    toUpdate = {}
    toAdd = UGDict.keys()
    for i in results:
        UniGeneID = i[0]
        if UGDict.has_key(UniGeneID):
            # If the key is in both the parsed file and the expression table, then it does not need to be added.
            toAdd.remove(i[0])
            # Compare the row from the expression table with the parsed file, to determine if the expression table needs updating.
            for j in columnIndices:
                if i[j] != UGDict[UniGeneID][j]:
                    if not toUpdate.has_key(UniGeneID):
                        toUpdate[UniGeneID] = []
                    toUpdate[UniGeneID].append(j)
        else:
            # If the key is in the expression table, but not in the parsed file, then the row needs to be removed.
            toRemove.append(i[0])
    values = '(' + ('%s,' * len(UGData[0]))
    values = values[:-1] + ')'
    
    #===========================================================================
    # Remove rows from the expression table that are not in the parsed file.
    #===========================================================================
    for i in toRemove:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.rowDELETE(cursor, tableUniGene, 'UniGeneID="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries removed from the UniGene table: ', len(toRemove)
    
    #===========================================================================
    # Update rows that have different values in the parsed file and the expression table.
    #===========================================================================
    for i in toUpdate.keys():
        toSet = []
        for j in toUpdate[i]:
            updateString = columns[j] + ' = "' + UGDict[i][j] + '"'
            toSet.append(updateString)
        toSet = ', '.join(toSet)
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableUPDATE(cursor, tableUniGene, toSet, 'UniGeneID="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries updated in the UniGene table: ', len(toUpdate)
    
    #===========================================================================
    # Add rows which are not in the expression table, but are in the parsed file.
    #===========================================================================
    rowsToAdd = [UGDict[i] for i in toAdd]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableINSERT(cursor, tableUniGene, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the UniGene table: ', len(toAdd)
    
    #===========================================================================
    # Enter the expression totals in the totals table.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('TRUNCATE TABLE ' + tableUniGeneTotals)
    values = '(' + ('%s,' * len(UGTotalsData[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableUniGeneTotals, values, UGTotalsData)
    mysql.closeConnection(conn, cursor)
    print '\tUniGene totals table updated.'