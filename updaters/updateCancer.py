import utilities.MySQLaccess as mysql

def main(CGCParsed, cancerTargets, UPExternalLinks, UPHumanAccessionMap, TTDTarget2Drug, DBTargetIDs, schemaProteins, tableCancerGene, DATABASEPASSWORD):

    # Record the mapping of non-representative UniProt accessions to representative accessions.
    accMap = {}
    readIn = open(UPHumanAccessionMap, 'r')
    for line in readIn:
        chunks = (line.strip()).split()
        accMap[chunks[0]] = chunks[1]
    readIn.close()

    # Read in the CGC information.
    CGCData = {}
    readIn = open(CGCParsed, 'r')
    for line in readIn:
        chunks = (line.strip()).split('\t')
        geneID = chunks[0]
        germline = chunks[1]
        somatic = chunks[2]
        CGCData[geneID] = {'Cancer' : 'Y', 'Target' : 'N', 'Somatic' : somatic, 'Germline' : germline}
    readIn.close()

    # Read in the information about anti-neoplastic drugs.
    cancerTTDTargets = set([])
    cancerDrugBankDrugs = set([])
    cancerDrugTargetUPAccs = set([])
    readIn = open(cancerTargets, 'r')
    readIn.readline()  # Strip the header line.
    for line in readIn:
        chunks = line.split('\t')
        cancerTTDTargets.add(chunks[2])
        cancerDrugBankDrugs.add(chunks[4])
    readIn.close()
    cancerTTDTargets -= set([''])
    cancerDrugBankDrugs -= set([''])

    # Determine the UniProt accessions of the TTD targets.
    readIn = open(TTDTarget2Drug, 'r')
    for line in readIn:
        accessions = line.split('\t')[1]
        for i in accessions.split(','):
            cancerDrugTargetUPAccs.add(i)
    readIn.close()

    # Determine the UniProt accessions of the targets of the DrugBank drugs.
    readIn = open(DBTargetIDs, 'r')
    for line in readIn:
        chunks = line.split('\t')
        accession = chunks[0]
        drugs = chunks[1].split(';')
        for i in drugs:
            if i in cancerDrugBankDrugs:
                cancerDrugTargetUPAccs.add(i)
    readIn.close()

    geneXRef = {}
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
    readIn.close()

    # Create the set of all representative UniProt accessions.
    allUPAccessions = set([])
    for i in accMap:
        allUPAccessions.add(accMap[i])

    # For every representative UniProt accession, determine whether it is implicated in cancer and if it is a target for anti-neoplastic drugs.
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
        # Record all the proteins that are targets of anti-neoplastic drugs as being implicated in cancer and being targets.
        try:
            # Only adding representative accessions.
            representativeAcc = accMap[i]
            cancerProteins[representativeAcc]['Cancer'] = 'Y'
            cancerProteins[representativeAcc]['Target'] = 'Y'
        except:
            pass

    dataRecords = [tuple([i, cancerProteins[i]['Cancer'], cancerProteins[i]['Target'], cancerProteins[i]['Somatic'], cancerProteins[i]['Germline']])
                   for i in cancerProteins]
    conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
    cursor.execute('TRUNCATE TABLE ' + tableCancerGene)
    values = '(' + ('%s,' * len(dataRecords[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableCancerGene, values, dataRecords)
    mysql.closeConnection(conn, cursor)