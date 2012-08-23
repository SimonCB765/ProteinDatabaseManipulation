'''
Created on 25 Oct 2011

@author: Simon Bull
'''

import utilities.MySQLaccess as mysql

def main(UPExternalLinks, tableUniProt2Ensembl, tableUniProt2GO, tableUniProt2UniGene, tableUniProt2HGNC,
         tableDict, schemaProteins, databasePassword):
    
    proteinLinks = {}
    readIn = open(UPExternalLinks, 'r')
    for line in readIn:
        line = line.strip()
        chunks = line.split(',')
        entrezLinks = [i for i in chunks[1].split(';')] if chunks[1] != '' else []
        unigeneLinks = [int(i.replace('Hs.', '')) for i in chunks[2].split(';')] if chunks[2] != '' else []
        goLinks = [int(i.replace('GO:', '')) for i in chunks[3].split(';')] if chunks[3] != '' else []
        hgncLinks = [i for i in chunks[4].split(';')] if chunks[4] != '' else []
        ensemblLinks = [i.split('-') for i in chunks[5].split(';')] if chunks[4] != '' else []
        proteinLinks[chunks[0]] = {'Entrez' : entrezLinks, 'UniGene' : unigeneLinks, 'GO' : goLinks,
                                   'HGNC' : hgncLinks, 'Ensembl' : ensemblLinks
                                   }
    readIn.close()
    
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    
    validXrefIDs = {'EnsemblGeneID' : set([]), 'GeneID' : set([]), 'GOTermID' : set([]), 'UPAccession' : set([]),
                    'UniGeneID' : set([]), 'EnsemblTranscriptID' : set([])}
    for i in tableDict.keys():
        for j in tableDict[i]:
            cursor = mysql.tableSELECT(cursor, i, j)
            results = [k[0] for k in cursor.fetchall()]
            results = set(results)
            validXrefIDs[i] = validXrefIDs[i].union(results)
    
    entrezInsert = []
    unigeneInsert = []
    goInsert = []
    hgncInsert = []
    ensemblInsert = []
    for i in proteinLinks.keys():
        if not i in validXrefIDs['UPAccession']:
            continue
        for j in proteinLinks[i]['Entrez']:
            if j in validXrefIDs['GeneID']:
                entrezInsert.append(tuple([i, str(j)]))
        for j in proteinLinks[i]['UniGene']:
            if j in validXrefIDs['UniGeneID']:
                unigeneInsert.append(tuple([i, j]))
        for j in proteinLinks[i]['GO']:
            if j in validXrefIDs['GOTermID']:
                goInsert.append(tuple([i, j]))
        for j in proteinLinks[i]['HGNC']:
            hgncInsert.append(tuple([i, j]))
        for j in proteinLinks[i]['Ensembl']:
            if j[0] in validXrefIDs['EnsemblGeneID']:# and j[1] in validXrefIDs['EnsemblTranscriptID']:
                ensemblInsert.append(tuple([i] + j))

##    print '\tNow Recording Entrez Crossreferences.'    
##    cursor.execute('TRUNCATE TABLE ' + tableUniProt2Entrez)
##    values = '(' + ('%s,' * len(entrezInsert[0]))
##    values = values[:-1] + ')'
##    mysql.tableINSERT(cursor, tableUniProt2Entrez, values, entrezInsert)

    print '\tNow Recording UniGene Crossreferences.'
    cursor.execute('TRUNCATE TABLE ' + tableUniProt2UniGene)
    values = '(' + ('%s,' * len(unigeneInsert[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableUniProt2UniGene, values, unigeneInsert)

    print '\tNow Recording GO Crossreferences.'
    cursor.execute('TRUNCATE TABLE ' + tableUniProt2GO)
    values = '(' + ('%s,' * len(goInsert[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableUniProt2GO, values, goInsert)

    print '\tNow Recording HGNC Crossreferences.'
    cursor.execute('TRUNCATE TABLE ' + tableUniProt2HGNC)
    values = '(' + ('%s,' * len(hgncInsert[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableUniProt2HGNC, values, hgncInsert)

    print '\tNow Recording Ensembl Crossreferences.'
    cursor.execute('TRUNCATE TABLE ' + tableUniProt2Ensembl)
    values = '(' + ('%s,' * len(ensemblInsert[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableUniProt2Ensembl, values, ensemblInsert)
    
    mysql.closeConnection(conn, cursor)