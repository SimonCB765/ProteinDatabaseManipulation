'''
Created on 17 Oct 2011

@author: Simon Bull
'''

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(parsedGOOutput, schemaProteins, tableGOInfo, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed GO data.
    #===========================================================================
    GOData = utilities.file2list.main(parsedGOOutput)
    GOData = [eval(i) for i in GOData]
    GODict = dict([(i[0], i) for i in GOData])
    
    #===========================================================================
    # Extract the GO information recorded in the database.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableSELECT(cursor, '*', tableGOInfo)
    results = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    
    #===========================================================================
    # Compare the parsed data with the data recorded in the table.
    #===========================================================================
    columnIndices = range(1, len(GOData[0]))
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('SHOW COLUMNS FROM ' + tableGOInfo)
    columns = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    columns = [i[0] for i in columns]
    
    toRemove = []
    toUpdate = {}
    toAdd = GODict.keys()
    for i in results:
        geneID = i[0]
        if GODict.has_key(geneID):
            # If the key is in both the parsed file and the table, then it does not need to be added.
            toAdd.remove(i[0])
            # Compare the row from the table with the parsed file, to determine if the table needs updating.
            for j in columnIndices:
                if i[j] != GODict[geneID][j]:
                    if not toUpdate.has_key(geneID):
                        toUpdate[geneID] = []
                    toUpdate[geneID].append(j)
        else:
            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
            toRemove.append(i[0])
    values = '(' + ('%s,' * len(GOData[0]))
    values = values[:-1] + ')'
    
    #===========================================================================
    # Remove rows from the table that are not in the parsed file.
    #===========================================================================
    for i in toRemove:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.rowDELETE(cursor, tableGOInfo, 'GOTermID="' + str(i) + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries removed from the GO table: ', len(toRemove)
    
    #===========================================================================
    # Update rows that have different values in the parsed file and the table.
    #===========================================================================
    for i in toUpdate.keys():
        toSet = []
        for j in toUpdate[i]:
            updateString = columns[j] + ' = "' + GODict[i][j] + '"'
            toSet.append(updateString)
        toSet = ', '.join(toSet)
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableUPDATE(cursor, tableGOInfo, toSet, 'GOTermID="' + str(i) + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries updated in the GO table: ', len(toUpdate)
    
    #===========================================================================
    # Add rows which are not in the table, but are in the parsed file.
    #===========================================================================
    # Split the records to be inserted into smaller chunks so the database connection is not lost
    # In order to get this to work you need to increase the size of max_allowed_packet for the
    # MySQL server. This is because the size of some of the paths is very large.
    # I did this by altering the default my.ini file to contain the line:
    # max_allowed_packets=32M
    # and put this line under the [mysqld] section
    rowsToAdd = [GODict[i] for i in toAdd]
    length = len(rowsToAdd)
    itemsInSplit = 5
    numberOfSplits = length / itemsInSplit
    if length%itemsInSplit != 0:
        numberOfSplits += 1
    recordsToInsert = [ rowsToAdd[i * itemsInSplit : (i+1) * itemsInSplit] for i in range(numberOfSplits)]
    for i in recordsToInsert:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableINSERT(cursor, tableGOInfo, values, i)
        mysql.closeConnection(conn, cursor)
    print '\tEntries added to the GO table: ', len(toAdd)