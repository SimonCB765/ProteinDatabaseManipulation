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

def update_homolog_table(ensemblParsedHomology, schemaProteins, tableHomologs, databasePassword):

    #===========================================================================
    # Extract and format the parsed gene data.
    #===========================================================================
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