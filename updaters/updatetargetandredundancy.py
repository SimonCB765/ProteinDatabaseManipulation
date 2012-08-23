'''
Created on 22 Oct 2011

@author: Simon Bull
'''

import culling.adjlistcreation
import culling.Leafcull

import utilities.file2list
import utilities.list2file
import utilities.MySQLaccess as mysql

def main(DBDrugIDs, DBTargetIDs, TTDUPAccessions, ChEMBLUPAccessions, UPHumanAccessionMap,
         UPDrugIDs, pharmGKBDrugTargets, folderCulling, schemaProteins, tableProteinInfo, tableNonRedundant, tableBLASTResults,
         databasePassword, viewsDict):
    
    allTargetsChEMBL = xref_chembl_uniprot(ChEMBLUPAccessions, UPHumanAccessionMap)
    allTargetsDB = xref_drugbank_uniprot(DBTargetIDs, UPHumanAccessionMap)
    allTargetsTTD = xref_TTD_uniprot(TTDUPAccessions, UPHumanAccessionMap)
    allTargetsUP = xref_uniprot_drugbank(UPDrugIDs, DBDrugIDs)
    allTargetPharmGKB = utilities.file2list.main(pharmGKBDrugTargets)
    print '\tTargets of approved drugs recorded by PharmGKB: ', len(allTargetPharmGKB)
    allTargets = list(set(allTargetsChEMBL) | set(allTargetsDB) | set(allTargetsTTD) | set(allTargetsUP) | set(allTargetPharmGKB))
    print '\tTotal number of unique targets found: ', len(allTargets)
    
    # Extract mode of action and clear the target information.
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    cursor = mysql.tableSELECT(cursor, 'UPAccession, ModeOfAction', tableProteinInfo)
    resultsModeOfAction = cursor.fetchall()
    for i in resultsModeOfAction:
        upid = i[0]
        mysql.tableUPDATE(cursor, tableProteinInfo, 'Target="N"', 'UPAccession="' + upid + '"')
##    cursor = mysql.tableSELECT(cursor, 'UPAccession, Target', tableProteinInfo)
##    resultsTargets = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    
    gpcr = []
    kinases = []
    ionChannels = []
    proteases = []
    for i in resultsModeOfAction:
        if i[1] == 'G-protein coupled receptor':
            gpcr.append(i[0])
        elif i[1] == 'Kinase':
            kinases.append(i[0])
        elif i[1] == 'Ion Channel':
            ionChannels.append(i[0])
        elif i[1] == 'Protease':
            proteases.append(i[0])

##    proteinTargets = {}
##    storedTargets = []
##    for i in resultsTargets:
##        target = False
##        if i[1] == 'Y':
##            target = True
##            storedTargets.append(i[0])
##        proteinTargets[i[0]] = target
##    allProteins = proteinTargets.keys()
    noLongerTargets = []
    newTargets = allTargets
##    storedNonTargets = [i for i in allProteins if i not in storedTargets]
##    noLongerTargets = list(set(storedTargets) - set(allTargets))
##    newTargets = list(set(allTargets) - set(storedTargets))
    
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    for i in noLongerTargets:
        mysql.tableUPDATE(cursor, tableProteinInfo, 'Target="N"', 'UPAccession="' + i + '"')
    for i in newTargets:
        mysql.tableUPDATE(cursor, tableProteinInfo, 'Target="Y"', 'UPAccession="' + i + '"')
    mysql.closeConnection(conn, cursor)
    
    # Check if the redundancy removal should be performed.
    if newTargets != [] or noLongerTargets != []:
        print '\tPerforming redundancy removal.'
        # Proteins have had their target status changed, so the redundancy needs to be recalculated.
        conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
        # Set the number of columns in the nonredundant table.
        cursor.execute('SHOW COLUMNS FROM ' + tableNonRedundant)
        numberColumns = len(cursor.fetchall())
        # Select all the proteins recorded in the database. The number of columns has one subtracted from it as the
        # UP accession column does not take the default 'NA' value.
        cursor = mysql.tableSELECT(cursor, 'UPAccession', tableProteinInfo)
        allProteins = [tuple([i[0]] + (['N'] * (numberColumns - 1))) for i in cursor.fetchall()]
        # Wipe and refill the nonredundant table.
        cursor.execute('TRUNCATE TABLE ' + tableNonRedundant)
        values = '(' + ('%s,' * numberColumns)
        values = values[:-1] + ')'
        mysql.tableINSERT(cursor, tableNonRedundant, values, allProteins)
        for column in sorted(viewsDict.keys()):
            print '\t\tRunning redundancy removal on ', column
            # For each set of proteins, run the Leaf program.
            inputLocation = folderCulling + '/' + column + '.txt'
            cursor = mysql.tableSELECT(cursor, '*', viewsDict[column])
            results = cursor.fetchall()
            # Determine the accessions of the proteins in the current set.
            proteinSet = [i[0] for i in results]
            print '\t\t\tSize of Redundant Dataset: ', len(proteinSet)
            proteinSetString = '\',\''.join(proteinSet)
            # Select all the BLAST results where both the hit and query protein are in the set to cull.
            cursor = mysql.tableSELECT(cursor, '*', tableBLASTResults, 'ProteinA IN (\'' + proteinSetString + '\') AND ProteinB IN (\'' + proteinSetString + '\')')
            protResults = cursor.fetchall()
            # Generate the file that is going to be used to perform the culling.
            writeTo = open(inputLocation, 'w')
            for i in protResults:
                writeTo.write('\t'.join([str(j) for j in i]) + '\n')
            writeTo.close()
            # Perform the culling.
            adjMatrix, proteinNames = culling.adjlistcreation.main(inputLocation, cutoffPercent=20, maxEValue=1, minAlignLength=20)
            print '\t\t\tNumber of Proteins in Similarity Graph: ', len(proteinNames)
            proteinsToCull = culling.Leafcull.main(adjMatrix, proteinNames)
            print '\t\t\tNumber of Proteins to Cull: ', len(proteinsToCull)
            for i in proteinsToCull:
                mysql.tableUPDATE(cursor, tableNonRedundant, column + '="N"', 'UPAccession="' + str(i) + '"')
            proteinsToKeep = [i for i in proteinSet if i not in proteinsToCull]
            print '\t\t\tNumber of Proteins to Keep: ', len(proteinsToKeep)
            for i in proteinsToKeep:
                mysql.tableUPDATE(cursor, tableNonRedundant, column + '="Y"', 'UPAccession="' + str(i) + '"')
        mysql.closeConnection(conn, cursor)

def xref_chembl_uniprot(ChEMBLUPAccessions, UPHumanAccessionMap):
    ChEMBlNonRepAccs = set([i.split('\t')[0] for i in utilities.file2list.main(ChEMBLUPAccessions)])
    UPAccMap = utilities.file2list.main(UPHumanAccessionMap)
    
    # Make a dictionary where the index is the, possibly deprecated UniProt accession, and the entry is the
    # representative accession.
    allUPID = {}
    for i in UPAccMap:
        chunks = i.split()
        allUPID[chunks[0]] = chunks[1]

    # Convert all the UniProt accessions recorded by ChEMBL, into the representative accessions.
    validUPLinks = [allUPID[i] for i in ChEMBlNonRepAccs if allUPID.has_key(i)]
    validUPLinks = list(set(validUPLinks))
    validUPLinks.sort()
    
    print '\tTargets of approved drugs recorded by ChEMBL: ', len(validUPLinks)
    
    return validUPLinks

def xref_drugbank_uniprot(DBTargetIDs, UPHumanAccessionMap):
    DBTargets = utilities.file2list.main(DBTargetIDs)
    UPAccMap = utilities.file2list.main(UPHumanAccessionMap)
    
    # Make a dictionary where the index is the, possibly deprecated UniProt accession, and the entry is the
    # representative accession.
    allUPID = {}
    for i in UPAccMap:
        chunks = i.split()
        allUPID[chunks[0]] = chunks[1]
    
    # Extract the UniProt accessions of the DrugBank targets.
    targetWithUPLink = [i.split('\t')[0] for i in DBTargets]
    # Determine which of the UniProt accessions from the DrugBank targets are for human proteins, and then convert all UniProt IDs
    # into the representative IDs.
    validUPLinks = [allUPID[i] for i in targetWithUPLink if allUPID.has_key(i)]
    validUPLinks = list(set(validUPLinks))
    validUPLinks.sort()
    
    print '\tTargets of approved drugs recorded by DrugBank: ', len(validUPLinks)
    
    return validUPLinks

def xref_TTD_uniprot(TTDUPAccessions, UPHumanAccessionMap):
    TTDNonRepAccs = utilities.file2list.main(TTDUPAccessions)
    UPAccMap = utilities.file2list.main(UPHumanAccessionMap)
    
    # Make a dictionary where the index is the, possibly deprecated UniProt accession, and the entry is the
    # representative accession.
    allUPID = {}
    for i in UPAccMap:
        chunks = i.split()
        allUPID[chunks[0]] = chunks[1]

    # Convert all the UniProt accessions recorded by the TTD, into the representative accessions.
    validUPLinks = [allUPID[i] for i in TTDNonRepAccs if allUPID.has_key(i)]
    validUPLinks = list(set(validUPLinks))
    validUPLinks.sort()
    
    print '\tTargets of approved drugs recorded by the TTD: ', len(validUPLinks)
    
    return validUPLinks

def xref_uniprot_drugbank(UPDrugIDs, DBDrugIDs):
    accsDrugsUP = utilities.file2list.main(UPDrugIDs)
    drugIDsDB = [(i.split('\t')) for i in utilities.file2list.main(DBDrugIDs)]
    drugIDsDB = set([i[0] for i in drugIDsDB if 'approved' in i[2]])
    
    approvedTargets = []
    for up in accsDrugsUP:
        chunks = up.split('\t')
        # If the DrugBank drug ID recorded by UniProt is in the list of approved DrugBank drug IDs,
        # record the UniProt accession as being a valid approved target.
        drugsUP = chunks[1].split(';')
        if drugIDsDB.intersection(drugsUP) != set([]):
            approvedTargets.append(chunks[0])

    # Ensure only unique results are returned.
    approvedTargets = list(set(approvedTargets))
    approvedTargets = [i for i in approvedTargets]
    approvedTargets.sort()

    print '\tTargets of approved drugs recorded by UniProt: ', str(len(approvedTargets))
    
    return approvedTargets