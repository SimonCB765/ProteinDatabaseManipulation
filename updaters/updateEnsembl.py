'''
Created on 18 Oct 2011

@author: Simon Bull
'''

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(ensemblParsedTranscripts, ensemblParsedGermVariants, ensemblParsedHomology,
         schemaProteins, tableEnsemblGene, tableGermVariants, tableHomologs, databasePassword):
    
    update_gene_table(ensemblParsedTranscripts, schemaProteins, tableEnsemblGene, databasePassword)
    update_variant_table(ensemblParsedGermVariants, schemaProteins, tableGermVariants, databasePassword, 'germ')
##    update_variant_table(ensemblParsedSomaticVariants, schemaProteins, tableCosmicVariants, databasePassword, 'cosmic')
    update_homolog_table(ensemblParsedHomology, schemaProteins, tableHomologs, databasePassword)

def update_gene_table(ensemblParsedTranscripts, schemaProteins, tableEnsemblGene, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed gene data.
    #===========================================================================
    geneData = utilities.file2list.main(ensemblParsedTranscripts)
    geneData = [eval(i) for i in geneData]
    geneDict = dict([(i[0], i) for i in geneData])
    values = '(' + ('%s,' * len(geneData[0]))
    values = values[:-1] + ')'
    
    rowsToAdd = geneData
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('TRUNCATE ' + tableEnsemblGene)
    cursor = mysql.tableINSERT(cursor, tableEnsemblGene, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the Ensembl gene table: ', len(rowsToAdd)

def update_variant_table(ensemblParsedGermVariants, schemaProteins, tableGermVariants, databasePassword, tableName):
    
    #===========================================================================
    # Extract and format the parsed variant data.
    #===========================================================================
##    variantData = utilities.file2list.main(ensemblParsedGermVariants)
##    variantData = [tuple([eval(j) if j.isdigit() else j for j in eval(i)]) for i in variantData]
##    variantData = [eval(i) for i in variantData]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('TRUNCATE ' + tableGermVariants)
    readVariants = open(ensemblParsedGermVariants, 'r')
    count = 0
    for line in readVariants:
        count += 1
        line = eval(line)
        values = '(' + ('%s,' * len(line))
        values = values[:-1] + ')'
        try:
            cursor = mysql.tableINSERT(cursor, tableGermVariants, values, [line])
        except:
            print line
            print count
            raise
    readVariants.close()
    mysql.closeConnection(conn, cursor)
##    variantDict = dict([((i[0], i[1]), i) for i in variantData])
##    
##    #===========================================================================
##    # Extract the variant information recorded in the database.
##    #===========================================================================
##    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##    cursor = mysql.tableSELECT(cursor, '*', tableGermVariants)
##    results = cursor.fetchall()
##    mysql.closeConnection(conn, cursor)
##    
##    #===========================================================================
##    # Compare the parsed data with the data recorded in the table.
##    #===========================================================================
##    columnIndices = range(1, len(variantData[0]))
##    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##    cursor.execute('SHOW COLUMNS FROM ' + tableGermVariants)
##    columns = cursor.fetchall()
##    mysql.closeConnection(conn, cursor)
##    columns = [i[0] for i in columns]
##    
##    toRemove = []
##    toUpdate = {}
##    toAdd = variantDict.keys()
##    for i in results:
##        transcriptID = i[0]
##        variantID = i[1]
##        primaryKey = (transcriptID, variantID)
##        if variantDict.has_key(primaryKey):
##            # If the key is in both the parsed file and the table, then it does not need to be added.
##            toAdd.remove(primaryKey)
##            # Compare the row from the table with the parsed file, to determine if the table needs updating.
##            for j in columnIndices:
##                if i[j] != variantDict[primaryKey][j]:
##                    if not toUpdate.has_key(primaryKey):
##                        toUpdate[primaryKey] = []
##                    toUpdate[primaryKey].append(j)
##        else:
##            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
##            toRemove.append(primaryKey)
##    values = '(' + ('%s,' * len(variantData[0]))
##    values = values[:-1] + ')'
##    
##    #===========================================================================
##    # Remove rows from the table that are not in the parsed file.
##    #===========================================================================
##    for i in toRemove:
##        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##        cursor = mysql.rowDELETE(cursor, tableGermVariants, 'EnsemblTranscriptID="' + i[0] + '" AND VariantID="' + i[1] + '"')
##        mysql.closeConnection(conn, cursor)
##    print '\tEntries removed from the ' + tableName + ' variant table: ', len(toRemove)
##    
##    #===========================================================================
##    # Update rows that have different values in the parsed file and the table.
##    #===========================================================================
##    for i in toUpdate.keys():
##        toSet = []
##        for j in toUpdate[i]:
##            updateString = columns[j] + ' = "' + str(variantDict[i][j]) + '"'
##            toSet.append(updateString)
##        toSet = ', '.join(toSet)
##        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##        cursor = mysql.tableUPDATE(cursor, tableGermVariants, toSet, 'EnsemblTranscriptID="' + i[0] + '" AND VariantID="' + i[1] + '"')
##        mysql.closeConnection(conn, cursor)
##    print '\tEntries updated in the ' + tableName + ' variant table: ', len(toUpdate)
##    
##    #===========================================================================
##    # Add rows which are not in the table, but are in the parsed file.
##    #===========================================================================
##    rowsToAdd = [variantDict[i] for i in toAdd]
##    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##    cursor = mysql.tableINSERT(cursor, tableGermVariants, values, rowsToAdd)
##    mysql.closeConnection(conn, cursor)
##    print '\tEntries added to the ' + tableName + ' variant table: ', len(toAdd)

def update_homolog_table(ensemblParsedHomology, schemaProteins, tableHomologs, databasePassword):
    
    #===========================================================================
    # Extract and format the parsed gene data.
    #===========================================================================
##    homologData = utilities.file2list.main(ensemblParsedHomology)
##    homologData = [eval(i) for i in homologData]
##    homologDict = dict([(tuple([i[0], i[1]]), i) for i in homologData])
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('TRUNCATE ' + tableHomologs)
    readVariants = open(ensemblParsedHomology, 'r')
    count = 0
    for line in readVariants:
        count += 1
        line = eval(line)
        values = '(' + ('%s,' * len(line))
        values = values[:-1] + ')'
        try:
            cursor = mysql.tableINSERT(cursor, tableHomologs, values, [line])
        except:
            print line
            print count
            raise
    readVariants.close()
    mysql.closeConnection(conn, cursor)
    
##    #===========================================================================
##    # Extract the gene information recorded in the database.
##    #===========================================================================
##    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##    cursor = mysql.tableSELECT(cursor, '*', tableHomologs)
##    results = cursor.fetchall()
##    mysql.closeConnection(conn, cursor)
##    
##    #===========================================================================
##    # Compare the parsed data with the data recorded in the table.
##    #===========================================================================
##    columnIndices = range(1, len(homologData[0]))
##    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##    cursor.execute('SHOW COLUMNS FROM ' + tableHomologs)
##    columns = cursor.fetchall()
##    mysql.closeConnection(conn, cursor)
##    columns = [i[0] for i in columns]
##    
##    toRemove = []
##    toUpdate = {}
##    toAdd = homologDict.keys()
##    for i in results:
##        humanGeneID = i[0]
##        homologGeneID = i[1]
##        dictKey = tuple([humanGeneID, homologGeneID])
##        if homologDict.has_key(dictKey):
##            # If the key is in both the parsed file and the table, then it does not need to be added.
##            toAdd.remove(dictKey)
##            # Compare the row from the table with the parsed file, to determine if the table needs updating.
##            for j in columnIndices:
##                if i[j] != homologDict[dictKey][j]:
##                    if not toUpdate.has_key(dictKey):
##                        toUpdate[dictKey] = []
##                    toUpdate[dictKey].append(j)
##        else:
##            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
##            toRemove.append(dictKey)
##    values = '(' + ('%s,' * len(homologData[0]))
##    values = values[:-1] + ')'
##    
##    #===========================================================================
##    # Remove rows from the table that are not in the parsed file.
##    #===========================================================================
##    for i in toRemove:
##        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##        cursor = mysql.rowDELETE(cursor, tableHomologs, 'HumanGene="' + i[0] + '" AND HomologGene="' + i[1] + '"')
##        mysql.closeConnection(conn, cursor)
##    print '\tEntries removed from the Ensembl homolog table: ', len(toRemove)
##    
##    #===========================================================================
##    # Update rows that have different values in the parsed file and the table.
##    #===========================================================================
##    for i in toUpdate.keys():
##        toSet = []
##        for j in toUpdate[i]:
##            updateString = columns[j] + ' = "' + str(homologDict[i][j]) + '"'
##            toSet.append(updateString)
##        toSet = ', '.join(toSet)
##        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##        cursor = mysql.tableUPDATE(cursor, tableHomologs, toSet, 'HumanGene="' + i[0] + '" AND HomologGene="' + i[1] + '"')
##        mysql.closeConnection(conn, cursor)
##    print '\tEntries updated in the Ensembl homolog table: ', len(toUpdate)
##    
##    #===========================================================================
##    # Add rows which are not in the table, but are in the parsed file.
##    #===========================================================================
##    rowsToAdd = [homologDict[i] for i in toAdd]
##    length = len(rowsToAdd)
##    itemsInSplit = 500
##    numberOfSplits = length / itemsInSplit
##    if length%itemsInSplit != 0:
##        numberOfSplits += 1
##    recordsToInsert = [ rowsToAdd[i * itemsInSplit : (i+1) * itemsInSplit] for i in range(numberOfSplits)]
##    count = 0
##    for i in recordsToInsert:
##        print count, numberOfSplits
##        count += 1
##        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
##        cursor = mysql.tableINSERT(cursor, tableHomologs, values, i)
##        mysql.closeConnection(conn, cursor)
##    print '\tEntries added to the Ensembl homolog table: ', len(toAdd)