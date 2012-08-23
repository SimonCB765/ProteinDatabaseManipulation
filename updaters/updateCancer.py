import utilities.MySQLaccess as mysql

def main(CGCParsed, cancerTargets, CGIHGNCIDs, CGIUPAccessions, UPExternalLinks, UPHumanAccessionMap,
         schemaProteins, tableCancerGene, tableCOSMICGene, tableCOSMICGene2Mutation, tableCOSMICMutation, DATABASEPASSWORD):

    accMap = {}
    readIn = open(UPHumanAccessionMap, 'r')
    for line in readIn:
        chunks = (line.strip()).split()
        accMap[chunks[0]] = chunks[1]
    readIn.close()

    CGCData = {}
    readIn = open(CGCParsed, 'r')
    for line in readIn:
        chunks = (line.strip()).split('\t')
        geneID = chunks[0]
        germline = chunks[1]
        somatic = chunks[2]
        CGCData[geneID] = {'Cancer' : 'Y', 'Target' : 'N', 'Somatic' : somatic, 'Germline' : germline}
    readIn.close()

    cancerDrugTargetUPAccs = set([])
    readIn = open(cancerTargets, 'r')
    readIn.readline()  # Strip the header line.
    for line in readIn:
        line = line[:-1]  # Strip the newline.
        UPAccs = line.split('\t')[10]
        UPAccs = set(UPAccs.split(';'))
        cancerDrugTargetUPAccs |= UPAccs
    readIn.close()
    cancerDrugTargetUPAccs -= set([''])

    CGIHGNCs = {}
    readIn = open(CGIHGNCIDs, 'r')
    for line in readIn:
        line = line.strip()
        CGIHGNCs[line] = {'Cancer' : 'Y', 'Target' : 'N', 'Somatic' : 'U', 'Germline' : 'U'}
    readIn.close()

    CGIUPAccs = {}
    readIn = open(CGIUPAccessions, 'r')
    for line in readIn:
        line = line.strip()
        CGIUPAccs[line] = {'Cancer' : 'Y', 'Target' : 'N', 'Somatic' : 'U', 'Germline' : 'U'}
    readIn.close()

    geneXRef = {}
    HGNCXRef = {}
    ensemblXRef = {}
    readIn = open(UPExternalLinks, 'r')
    for line in readIn:
        chunks = line.split(',')
        UPAcc = chunks[0]
        geneID = chunks[1]
        if geneID == '':
            pass
        else:
            geneID = geneID.split(';')
            for i in geneID:
                if geneXRef.has_key(i):
                    geneXRef[i].add(UPAcc)
                else:
                    geneXRef[i] = set([UPAcc])
        hgnc = chunks[4]
        if hgnc == '':
            pass
        else:
            hgnc = hgnc.split(';')
            for i in hgnc:
                if HGNCXRef.has_key(i):
                    HGNCXRef[i].add(UPAcc)
                else:
                    HGNCXRef[i] = set([UPAcc])
        ensemblTrans = chunks[5]
        if ensemblTrans == '':
            pass
        else:
            ensemblTrans = [i.split('-')[1] for i in ensemblTrans.split(';')]
            for i in ensemblTrans:
                if ensemblXRef.has_key(i):
                    ensemblXRef[i].add(UPAcc)
                else:
                    ensemblXRef[i] = set([UPAcc])
    readIn.close()

##    conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
##    # Get all COSMIC genes that have a mutation.
##    cursor = mysql.tableSELECT(cursor, '*', tableCOSMICGene2Mutation)
##    results = cursor.fetchall()
##    mutatedGenes = set([i[0] for i in results])
##    mutation2gene = {}
##    for i in results:
##        gene = i[0]
##        mutation = i[1]
##        if mutation2gene.has_key(mutation):
##            mutation2gene[mutation].add(gene)
##        else:
##            mutation2gene[mutation] = set([gene])
##
##    # Get the HGNC and Ensembl Transcript IDs for the mutated genes.
##    cursor = mysql.tableSELECT(cursor, 'gene, transcript, HGNC', tableCOSMICGene, 'gene IN (\'' + '\', \''.join(mutatedGenes) + '\')')
##    results = cursor.fetchall()
##    gene2IDs = dict([(i[0], [i[1], i[2]]) for i in results])
##
##    # Get the information about the type of mutations.
##    cursor = mysql.tableSELECT(cursor, 'mutation, mutationtype', tableCOSMICMutation, 'mutation in (\'' + '\', \''.join([str(i) for i in mutation2gene.keys()]) + '\')')
##    results = cursor.fetchall()
##    mut2type = dict([(i[0], i[1]) for i in results])
##
##    mysql.closeConnection(conn, cursor)

##    # COSMIC genes invloved in cancer are taken to be genes that are mutated in a cancerous sample.
##    gene2type = {}
##    for i in mutation2gene.keys():
##        gene = mutation2gene[i]
##        mutType = mut2type[i]
##        for i in gene:
##            if gene2type.has_key(i):
##                gene2type[i][mutType] = True
##            else:
##                gene2type[i] = {'Somatic' : False, 'Germline' : False, 'Unknown' : False}
##                gene2type[i][mutType] = True
##
##    ensemblTransCancerTypes = {}
##    HGNCCancerTypes = {}
##    for i in gene2type.keys():
##        ensemblTrans = gene2IDs[i][0]
##        if ensemblTrans != '':
##            if ensemblTransCancerTypes.has_key(ensemblTrans):
##                for j in gene2type[i]:
##                    ensemblTransCancerTypes[ensemblTrans][j] = ensemblTransCancerTypes[ensemblTrans][j] or gene2type[i][j]
##            else:
##                ensemblTransCancerTypes[ensemblTrans] = gene2type[i]
##
##        hgnc = gene2IDs[i][1]
##        if hgnc != -1:
##            hgnc = str(hgnc)
##            if HGNCCancerTypes.has_key(hgnc):
##                for j in gene2type[i]:
##                    HGNCCancerTypes[hgnc][j] = HGNCCancerTypes[hgnc][j] or gene2type[i][j]
##            else:
##                HGNCCancerTypes[hgnc] = gene2type[i]

    allUPAccessions = set([])
    for i in accMap:
        allUPAccessions.add(accMap[i])
    cancerProteins = dict([(i, {'Cancer' : 'N', 'Target' : 'N', 'Somatic' : 'N', 'Germline' : 'N'}) for i in allUPAccessions])
    for i in CGCData.keys():
        if geneXRef.has_key(i):
            for j in geneXRef[i]:
                cancerProteins[j]['Cancer'] = 'Y'
                if cancerProteins[j]['Somatic'] == 'Y':
                    # Keep it as a yes.
                    pass
                elif cancerProteins[j]['Somatic'] == 'N':
                    if CGCData[i]['Somatic'] == 'Y':
                        # Change it to a yes.
                        cancerProteins[j]['Somatic'] = 'Y'
                if cancerProteins[j]['Germline'] == 'Y':
                    # Keep it as a yes.
                    pass
                elif cancerProteins[j]['Germline'] == 'N':
                    if CGCData[i]['Germline'] == 'Y':
                        # Change it to a yes.
                        cancerProteins[j]['Germline'] = 'Y'

    for i in cancerDrugTargetUPAccs:
        cancerProteins[i]['Cancer'] = 'Y'
        cancerProteins[i]['Target'] = 'Y'

##    for i in CGIHGNCs.keys():
##        if HGNCXRef.has_key(i):
##            for j in HGNCXRef[i]:
##                cancerProteins[j]['Cancer'] = 'Y'
##                if cancerProteins[j]['Somatic'] == 'Y':
##                    # Keep it as a yes.
##                    pass
##                else:
##                    cancerProteins[j]['Somatic'] = 'U'
##                if cancerProteins[j]['Germline'] == 'Y':
##                    # Keep it as a yes.
##                    pass
##                else:
##                    cancerProteins[j]['Germline'] = 'U'

##    for i in CGIUPAccs.keys():
##        if accMap.has_key(i):
##            reprUPAcc = accMap[i]
##            cancerProteins[reprUPAcc]['Cancer'] = 'Y'
##            if cancerProteins[j]['Somatic'] == 'Y':
##                # Keep it as a yes.
##                pass
##            else:
##                cancerProteins[j]['Somatic'] = 'U'
##            if cancerProteins[j]['Germline'] == 'Y':
##                # Keep it as a yes.
##                pass
##            else:
##                cancerProteins[j]['Germline'] = 'U'

##    for i in ensemblTransCancerTypes.keys():
##        if ensemblXRef.has_key(i):
##            for j in ensemblXRef[i]:
##                cancerProteins[j]['Cancer'] = 'Y'
##                if cancerProteins[j]['Somatic'] == 'Y' or ensemblTransCancerTypes[i]['Somatic']:
##                    cancerProteins[j]['Somatic'] = 'Y'
##                elif ensemblTransCancerTypes[i]['Unknown']:
##                    cancerProteins[j]['Somatic'] = 'U'
##                if cancerProteins[j]['Germline'] == 'Y' or ensemblTransCancerTypes[i]['Germline']:
##                    cancerProteins[j]['Germline'] = 'Y'
##                elif ensemblTransCancerTypes[i]['Unknown']:
##                    cancerProteins[j]['Germline'] = 'U'
##
##    for i in HGNCCancerTypes.keys():
##        if HGNCXRef.has_key(i):
##            for j in HGNCXRef[i]:
##                cancerProteins[j]['Cancer'] = 'Y'
##                if cancerProteins[j]['Somatic'] == 'Y' or HGNCCancerTypes[i]['Somatic']:
##                    cancerProteins[j]['Somatic'] = 'Y'
##                elif HGNCCancerTypes[i]['Unknown']:
##                    cancerProteins[j]['Somatic'] = 'U'
##                if cancerProteins[j]['Germline'] == 'Y' or HGNCCancerTypes[i]['Germline']:
##                    cancerProteins[j]['Germline'] = 'Y'
##                elif HGNCCancerTypes[i]['Unknown']:
##                    cancerProteins[j]['Germline'] = 'U'

    dataRecords = [tuple([i, cancerProteins[i]['Cancer'], cancerProteins[i]['Target'], cancerProteins[i]['Somatic'], cancerProteins[i]['Germline']])
                   for i in cancerProteins]
    conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
    cursor.execute('TRUNCATE TABLE ' + tableCancerGene)
    values = '(' + ('%s,' * len(dataRecords[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableCancerGene, values, dataRecords)
    mysql.closeConnection(conn, cursor)