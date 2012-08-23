'''
Created on 10 Nov 2011

@author: Simon Bull
'''

import utilities.MySQLaccess as mysql
import utilities.file2list

def main(UPPPIData, schemaProteins, tablePPI, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed gene data.
    #===========================================================================
    ppiData = utilities.file2list.main(UPPPIData)
    ppiData = [eval(i) for i in ppiData]
    ppiDict = dict([(tuple([i[0], i[1]]), i) for i in ppiData])
    
    #===========================================================================
    # Extract the gene information recorded in the database.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableSELECT(cursor, '*', tablePPI)
    results = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    
    #===========================================================================
    # Compare the parsed data with the data recorded in the table.
    #===========================================================================
    columnIndices = range(1, len(ppiData[0]))
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('SHOW COLUMNS FROM ' + tablePPI)
    columns = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    columns = [i[0] for i in columns]
    
    toRemove = []
    toUpdate = {}
    toAdd = ppiDict.keys()
    for i in results:
        proteinOne = i[0]
        proteinTwo = i[1]
        dictKey = tuple([proteinOne, proteinTwo])
        if ppiDict.has_key(dictKey):
            # If the key is in both the parsed file and the table, then it does not need to be added.
            toAdd.remove(dictKey)
            # Compare the row from the table with the parsed file, to determine if the table needs updating.
            for j in columnIndices:
                if i[j] != ppiDict[dictKey][j]:
                    if not toUpdate.has_key(dictKey):
                        toUpdate[dictKey] = []
                    toUpdate[dictKey].append(j)
        else:
            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
            toRemove.append(dictKey)
    values = '(' + ('%s,' * len(ppiData[0]))
    values = values[:-1] + ')'
    
    #===========================================================================
    # Remove rows from the table that are not in the parsed file.
    #===========================================================================
    for i in toRemove:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.rowDELETE(cursor, tablePPI, 'PPIProteinOne="' + i[0] + '" AND PPIProteinTwo="' + i[1] + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries removed from the PPI table: ', len(toRemove)
    
    #===========================================================================
    # Update rows that have different values in the parsed file and the table.
    #===========================================================================
    for i in toUpdate.keys():
        toSet = []
        for j in toUpdate[i]:
            updateString = columns[j] + ' = "' + str(ppiDict[i][j]) + '"'
            toSet.append(updateString)
        toSet = ', '.join(toSet)
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableUPDATE(cursor, tablePPI, toSet, 'PPIProteinOne="' + i[0] + '" AND PPIProteinTwo="' + i[1] + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries updated in the PPI table: ', len(toUpdate)
    
    #===========================================================================
    # Add rows which are not in the table, but are in the parsed file.
    #===========================================================================
    rowsToAdd = [ppiDict[i] for i in toAdd]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableINSERT(cursor, tablePPI, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the PPI table: ', len(toAdd)