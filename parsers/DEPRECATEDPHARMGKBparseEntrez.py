import re
import utilities.MySQLaccess as mysql

def main(entrezHumanGeneIDs, CGCData, MeSHData, UPExternalLinks, entrezParsedOutput, pharmGKBDrugTargets,
         schemaPharmGKB, databasePassword):

    parsedCGC = parse_CGC(CGCData)
    MeSHCancerTerms, MeSHCancerIDs = parse_MeSH(MeSHData, 'C04')  # C04 is the root term for neoplasms.
    MeSHCancerIDs = '(\'' + '\', \''.join(MeSHCancerIDs) + '\')'

    # Access the schema containing the PharmGKB data.
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaPharmGKB)
    cursor = mysql.tableSELECT(cursor, 'disease', schemaPharmGKB + '.disease2mesh', 'mesh IN ' + MeSHCancerIDs)
    pharmgkbCancerDiseaseIDs = [i[0] for i in cursor.fetchall()]  # PharmGKB IDs for disease linked to a neoplasm MeSH term.
    pharmgkbCancerDiseaseIDs = '(\'' + '\', \''.join(pharmgkbCancerDiseaseIDs) + '\')'

    cursor = mysql.tableSELECT(cursor, 'gene', schemaPharmGKB + '.gene2disease', 'disease IN ' + pharmgkbCancerDiseaseIDs)
    pharmgkkbCancerGenes = [i[0] for i in cursor.fetchall()]  # PharmGKB IDs for genes linked to a neoplasm MeSH term.

    cursor = mysql.tableSELECT(cursor, 'drug', schemaPharmGKB + '.drug2disease', 'disease IN ' + pharmgkbCancerDiseaseIDs)
    pharmgkkbCancerDrugs = [i[0] for i in cursor.fetchall()]  # PharmGKB IDs for drugs linked to a neoplasm MeSH term.
##    cursor = mysql.tableSELECT(cursor, 'drugclass', schemaPharmGKB + '.drugclass2disease', 'disease IN ' + pharmgkbCancerDiseaseIDs)
##    pharmgkkbCancerDrugClasses = [i[0] for i in cursor.fetchall()]  # PharmGKB IDs for drug classes linked to a neoplasm MeSH term.

    cursor = mysql.tableSELECT(cursor, 'gene', schemaPharmGKB + '.drug2gene', 'drug IN ' + '(\'' + '\', \''.join(pharmgkkbCancerDrugs) + '\')')
    pharmgkkbDrugTargetsOfCancerDrugs = set([i[0] for i in cursor.fetchall()])  # PharmGKB IDs for genes linked to a drug linked to a neoplasm MeSH term.
##    cursor = mysql.tableSELECT(cursor, 'gene', schemaPharmGKB + '.drugclass2gene', 'drugclass IN ' + '(\'' + '\', \''.join(pharmgkkbCancerDrugClasses) + '\')')
##    pharmgkkbDrugClassTargetsOfCancerDrugs = set([i[0] for i in cursor.fetchall()])  # PharmGKB IDs for genes linked to a drug class linked to a neoplasm MeSH term.

    pharmgkbCancerDrugTargets = pharmgkkbDrugTargetsOfCancerDrugs.intersection(pharmgkkbCancerGenes)  # PharmGKB IDs for genes linked to a neoplasm MeSH term and a drug linked to a neoplasm MeSH term.
##    pharmgkbCancerDrugClassTargets = pharmgkkbDrugClassTargetsOfCancerDrugs.intersection(pharmgkkbCancerGenes)  # PharmGKB IDs for genes linked to a neoplasm MeSH term and a drug class linked to a neoplasm MeSH term.

##    cursor = mysql.tableSELECT(cursor, 'entrezid', schemaPharmGKB + '.gene', 'geneid IN ' + '(\'' + '\', \''.join(pharmgkbCancerDrugTargets.union(pharmgkbCancerDrugClassTargets)) + '\')')
    cursor = mysql.tableSELECT(cursor, 'entrezid', schemaPharmGKB + '.gene', 'geneid IN ' + '(\'' + '\', \''.join(pharmgkbCancerDrugTargets) + '\')')
    entrezPharmGKBCancerDrugTargets = sorted([int(i[0]) for i in cursor.fetchall()])  # Entrez Gene IDs for targets of drugs linked to a neoplasm MeSH term.
    cursor = mysql.tableSELECT(cursor, 'entrezid', schemaPharmGKB + '.gene', 'geneid IN ' + '(\'' + '\', \''.join(pharmgkkbCancerGenes) + '\')')
    entrezPharmGKBCancerNonDrugTargets = set([i[0] for i in cursor.fetchall()])
    entrezPharmGKBCancerNonDrugTargets = entrezPharmGKBCancerNonDrugTargets.union(parsedCGC.keys())
    entrezPharmGKBCancerNonDrugTargets = sorted([int(i) for i in entrezPharmGKBCancerNonDrugTargets if not i in entrezPharmGKBCancerDrugTargets])  # Entrez Gene IDs for non-targets of drugs linked to a neoplasm MeSH term.
    mysql.closeConnection(conn, cursor)

    parsedHumanGenes = sorted(parse_human_genes(entrezHumanGeneIDs))

    writeTo = open(entrezParsedOutput, 'w')
    for i in parsedHumanGenes:
        # Line format is:
        # GeneID\tCancer\tTarget\tSomatic\tGermline\n
        outputLine = str(i)
        if i in parsedCGC.keys():
            if i in entrezPharmGKBCancerDrugTargets:
                outputLine += '\tY\tY'
            else:
                outputLine += '\tY\tN'
            outputLine += '\t' + parsedCGC[i]['Somatic'] + '\t' + parsedCGC[i]['Germline'] + '\n'
        else:
            if i in entrezPharmGKBCancerDrugTargets:
                outputLine += '\tY\tY\tU\tU\n'
            elif i in entrezPharmGKBCancerNonDrugTargets:
                outputLine += '\tY\tN\tU\tU\n'
            else:
                outputLine += '\tN\tN\tN\tN\n'
        writeTo.write(outputLine)
    writeTo.close()

    # Extract all drug targets from PharmGKB.
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaPharmGKB)
    cursor = mysql.tableSELECT(cursor, 'DISTINCT gene', schemaPharmGKB + '.drug2gene', 'pathwayevidence = \'Y\'')
    drugGenes = [i[0] for i in cursor.fetchall()]
##    cursor = mysql.tableSELECT(cursor, 'DISTINCT gene', schemaPharmGKB + '.drugclass2gene')
##    drugclassGenes = cursor.fetchall()
##    geneTargets = [i[0] for i in set(drugGenes).union(drugclassGenes)]
    cursor = mysql.tableSELECT(cursor, 'DISTINCT entrezid', schemaPharmGKB + '.gene', 'geneid IN (\'' + '\', \''.join(drugGenes) + '\')')
    entrezTargets = [i[0] for i in cursor.fetchall()]
    mysql.closeConnection(conn, cursor)

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

    UPDrugTargets = set([])
    for i in entrezTargets:
        i = str(i)
        if geneXRef.has_key(i):
            UPDrugTargets.update(geneXRef[i])

    writeTo = open(pharmGKBDrugTargets, 'w')
    for i in sorted(UPDrugTargets):
        pass#writeTo.write(i + '\n')
    writeTo.close()
    

def parse_CGC(CGCData):

    readIn = open(CGCData, 'r')
    data = readIn.read()
    readIn.close()
    data = data.split('\r')[1:]  # For some reason the file is split on carriage returns, so the normal file reading method doesn't work. [1:] in order to skip the header line.
    parsedCGC = {}
    for line in data:
        chunks = line.split('\t')
        ncbiGene = int(chunks[2])
        somatic = 'Y' if 'yes' in chunks[5] else 'N'
        germline = 'Y' if 'yes' in chunks[6] else 'N'
        if parsedCGC.has_key(ncbiGene):
            somatic = 'Y' if somatic == 'Y' or parsedCGC[ncbiGene]['Somatic'] == 'Y' else 'N'
            germline = 'Y' if germline == 'Y' or parsedCGC[ncbiGene]['Germline'] == 'Y' else 'N'
        parsedCGC[ncbiGene] = {'Somatic' : somatic, 'Germline' : germline}

    return parsedCGC

def parse_MeSH(MeSHData, rootTerm):
    neoplasmTerms = set([])
    neoplasmIDs = set([])

    readIn = open(MeSHData, 'r')
    validRecord = False
    term = ''
    ID = ''
    for line in readIn:
        line = line.rstrip()
        if line == '*NEWRECORD':
            if validRecord:
                neoplasmTerms.add(term)
                neoplasmIDs.add(ID)
            validRecord = False
            term = ''
            ID = ''
        elif line[:5] == 'MH = ':
            term = line.split(' = ')[1]
        elif line[:5] == 'MN = ':
            validRecord = True if re.search('^' + rootTerm, line.split(' = ')[1]) else False
        elif line[:5] == 'UI = ':
            ID = line.split(' = ')[1]
    readIn.close()

    neoplasmTerms -= set([''])
    neoplasmIDs -= set([''])

    return neoplasmTerms, neoplasmIDs

def parse_human_genes(entrezHumanGeneIDs):

    humanGenes = set([])
    readIn = open(entrezHumanGeneIDs, 'r')
    header = readIn.readline()
    for line in readIn:
        chunks = line.split('\t')
        humanGenes.add(int(chunks[1]))
    readIn.close()

    return humanGenes