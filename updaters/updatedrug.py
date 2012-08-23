'''
Created on 17 Nov 2011

@author: Simon Bull
'''

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(UPDrugIDs, DBDrugIDs, DBTargetIDs, TTDTarget2Drug, ChEMBLUPAccessions, ChEMBLCID, bindingParsed,
         UPHumanAccessionMap, databasePassword, schemaProteins, tableDrugs):
    
    # Generate the human accession mapping:
    UPAccMap = utilities.file2list.main(UPHumanAccessionMap)
    # Make a dictionary where the index is the, possibly deprecated UniProt accession, and the entry is the
    # representative accession.
    allUPID = {}
    for i in UPAccMap:
        chunks = i.split()
        allUPID[chunks[0]] = chunks[1]
    
    # Extract the DrugBank drug IDs, the CAS number and the CID for every approved DrugBank drug.
    approvedDrugBankDrugs = {}
    approvedDrugTuples = set([])
    readDrugs = open(DBDrugIDs, 'r')
    for line in readDrugs:
        chunks = (line.replace('\n', '')).split('\t')
        if len(chunks) < 4 or not 'approved' in chunks[2]:
            continue
        name = chunks[1]
        CAS = chunks[3]
        CID = chunks[4]
        approvedDrugBankDrugs[chunks[0]] = tuple([name, CAS, CID])
        approvedDrugTuples.add(tuple([name, CAS, CID]))
    readDrugs.close()
    approvedDrugs = approvedDrugBankDrugs.keys()
    
#    print len(approvedDrugs)
    
    # Extract the information about which DrugBank drugs target UniProt proteins, as recorded by UniProt.
    targetsUP = utilities.file2list.main(UPDrugIDs)
    targetDrugLinks = {}
    for i in targetsUP:
        chunks = i.split('\t')
        targetDrugLinks[chunks[0]] = set([i for i in chunks[1].split(';') if i in approvedDrugs])
    
    # Extract the information about which DrugBank drugs target UniProt proteins, as recorded by DrugBank.
    targetsDB = utilities.file2list.main(DBTargetIDs)
    for i in targetsDB:
        chunks = i.split('\t')
        proteinAcc = chunks[0]
        if not allUPID.has_key(proteinAcc):
            continue
        proteinAcc = allUPID[proteinAcc]
        if targetDrugLinks.has_key(proteinAcc):
            targetDrugLinks[proteinAcc].union([i for i in chunks[1].split(';') if i in approvedDrugs])
            #targetDrugLinks[proteinAcc].intersection([i for i in chunks[1].split(';') if i in approvedDrugs])
        else:
            targetDrugLinks[proteinAcc] = [i for i in chunks[1].split(';') if i in approvedDrugs]
    
    # Convert DrugBank drug IDs to CIDs where possible.
    for i in targetDrugLinks.keys():
        targetDrugLinks[i] = [approvedDrugBankDrugs[j] for j in targetDrugLinks[i]]
    
#    for i in targetDrugLinks.keys()[:50]:
#        print i, len(targetDrugLinks[i])
#    print len(targetDrugLinks.keys()), sorted(targetDrugLinks.keys())
    
    # Extract the TTD target-drug relationship information.
    targetsTTD = utilities.file2list.main(TTDTarget2Drug)
    for i in targetsTTD:
        chunks = i.split('\t')
        proteinAccs = [allUPID[i] for i in chunks[1].split(',') if allUPID.has_key(i)]
        if proteinAccs == []:
            continue
        drugInfo = [tuple(j.split(',')) for j in chunks[2].split(';')]
        for j in drugInfo:
            approvedDrugTuples.add(j)
        for j in proteinAccs:
            if not targetDrugLinks.has_key(j):
                # If the UniProt accession associated with the TTD target is not already recorded as a drug target,
                # then add it along with the drugs that target it.
                targetDrugLinks[j] = drugInfo
            else:
                # The UniProt accession is already recorded in the target drug relationship dictionary.
                # Collect all the (Name, CAS, CID) tuples, and see if the drugs associated with the UniProt
                # accession in the TTD, are already accounted for.
                drugTuples = targetDrugLinks[j]
                drugNames = set([])
                drugCASes = set([])
                drugCIDs = set([])
                for k in drugTuples:
                    drugNames.add(k[0])
                    drugCASes.add(k[1])
                    drugCIDs.add(k[2])
                emptySet = set([''])
                drugNames -= emptySet
                drugCASes -= emptySet
                drugCIDs -= emptySet
                for k in drugInfo:
                    if k[0] in drugNames or k[1] in drugCASes or k[2] in drugCIDs:
                    # If this is True, then that means that one of the name, CAS or CID is already recorded in the
                    # list of drugs that target the protein. Therefore, the drug is already recorded as targeting
                    # the protein and should not be added.
                        continue
                    targetDrugLinks[j].append(k)
    
#    for i in targetDrugLinks.keys()[:50]:
#        print i, len(targetDrugLinks[i])
#    print len(targetDrugLinks.keys()), sorted(targetDrugLinks.keys())
    
    # Extract the ChEMBL drug-target relationship information, along with the binding information.
    molregno2CID = {}
    readCIDs = open(ChEMBLCID, 'r')
    for line in readCIDs:
        chunks = (line.strip()).split('\t')
        molregno2CID[chunks[0]] = chunks[1]
    readCIDs.close()
    targetDrugBinding = {}
    targetsChEMBL = utilities.file2list.main(ChEMBLUPAccessions)
    for i in targetsChEMBL:
        chunks = i.split('\t')
        UPAccession = chunks[0]
        if not allUPID.has_key(UPAccession):
            # If the target protein is not a UniProt human protein.
            continue
        else:
            UPAccession = allUPID[UPAccession]
        molregno = chunks[1]
        drugID = molregno2CID[molregno] if molregno2CID.has_key(molregno) else 'm' + molregno
        name = chunks[2]
        
        # Update the drug-target relationship information.
        drugInfo = tuple([name, '', drugID])
        if not targetDrugLinks.has_key(UPAccession):
            # If the UniProt accession associated with the TTD target is not already recorded as a drug target,
            # then add it along with the drugs that target it.
            targetDrugLinks[UPAccession] = [drugInfo]
        else:
            # The UniProt accession is already recorded in the target drug relationship dictionary.
            # Collect all the (Name, CAS, CID) tuples, and see if the drugs associated with the UniProt
            # accession in the TTD, are already accounted for.
            drugTuples = targetDrugLinks[UPAccession]
            drugNames = set([])
            drugCASes = set([])
            drugCIDs = set([])
            for j in drugTuples:
                drugNames.add(j[0])
                drugCASes.add(j[1])
                drugCIDs.add(j[2])
            emptySet = set([''])
            drugNames -= emptySet
            drugCASes -= emptySet
            drugCIDs -= emptySet
            if drugInfo[0] in drugNames or drugInfo[1] in drugCASes or drugInfo[2] in drugCIDs:
            # If this is True, then that means that one of the name, CAS or CID is already recorded in the
            # list of drugs that target the protein. Therefore, the drug is already recorded as targeting
            # the protein and should not be added.
                pass
            else:
                targetDrugLinks[UPAccession].append(drugInfo)
        
        # Record drug binding information.
        if chunks[4] == 'None' or chunks[5] == 'None':
            continue
        value = float(chunks[4])
        units = chunks[5]
        if units == 'M':
            # Convert the value to nM.
            value *= 1000000000
        elif units == 'mM':
            # Convert the value to nM.
            value *= 1000000
        elif units == 'uM':
            # Convert the value to nM.
            value *= 1000
        type = chunks[6]
        bindingKey = tuple([UPAccession, drugID])
        if type.lower() == 'ki':
            if not targetDrugBinding.has_key(bindingKey):
                targetDrugBinding[bindingKey] = {'Ki' : value, 'Kd' : float('Inf')}
            else:
                targetDrugBinding[bindingKey]['Ki'] = min(targetDrugBinding[bindingKey]['Ki'], value)
        elif type.lower() == 'kd':
            if not targetDrugBinding.has_key(bindingKey):
                targetDrugBinding[bindingKey] = {'Ki' : float('Inf'), 'Kd' : value}
            else:
                targetDrugBinding[bindingKey]['Kd'] = min(targetDrugBinding[bindingKey]['Kd'], value)

#    for i in targetDrugLinks.keys()[:50]:
#        print i, len(targetDrugLinks[i]), targetDrugLinks[i]
#    print len(targetDrugLinks.keys()), sorted(targetDrugLinks.keys())
    
    # Extract binding information from BindingDB.
    bindingData = utilities.file2list.main(bindingParsed)
    for i in bindingData:
        chunks = i.split('\t')
        UPAccession = chunks[0]
        if allUPID.has_key(UPAccession):
            UPAccession = allUPID[UPAccession]
        else:
            continue
        CID = chunks[1]
        Ki = chunks[2]
        Kd = chunks[3]
        if CID == 'n/a' or (Ki == 'n/a' and Kd == 'n/a'):
            # If any of these are True then the binding information can not be cross-referenced (CID == 'n/a',
            # or there is no binding information.
            continue
        if Ki[0] in ['>', '<', '=']:
            Ki = float(Ki[1:])
        else:
            try:
                Ki = float(Ki)
            except:
                pass
        if Kd[0] in ['>', '<', '=']:
            Ki = float(Kd[1:])
        else:
            try:
                Kd = float(Kd)
            except:
                pass
        bindingKey = tuple([UPAccession, CID])
        if Ki != 'n/a':
            if not targetDrugBinding.has_key(bindingKey):
                targetDrugBinding[bindingKey] = {'Ki' : Ki, 'Kd' : float('Inf')}
            else:
                targetDrugBinding[bindingKey]['Ki'] = min(targetDrugBinding[bindingKey]['Ki'], Ki)
        if Kd != 'n/a':
            if not targetDrugBinding.has_key(bindingKey):
                targetDrugBinding[bindingKey] = {'Ki' : float('Inf'), 'Kd' : Kd}
            else:
                targetDrugBinding[bindingKey]['Kd'] = min(targetDrugBinding[bindingKey]['Kd'], Kd)
    
    # Generate the tuples to insert into the database.
    tuplesToInsert = []
    for i in targetDrugLinks:
        for j in targetDrugLinks[i]:
            drugName = j[0]
            drugID = j[2]
            bindingKey = tuple([i, drugID])
            if targetDrugBinding.has_key(bindingKey):
                bindingData = targetDrugBinding[bindingKey]
                Ki = bindingData['Ki']
                if Ki == float('inf'):
                    Ki = -1
                Kd = bindingData['Kd']
                if Kd == float('inf'):
                    Kd = -1
            else:
                Ki = -1
                Kd = -1
            tuplesToInsert.append(tuple([i, drugID, drugName, Ki, Kd]))
    
    # Remove potential duplicates caused by the names of the drugs having different capital letters, and mysql
    # ignoring the capitalisation of the names.
    tuplesToInsert = [tuple([i[0], i[1], i[2].upper(), i[3], i[4]]) for i in tuplesToInsert]
    tuplesToInsert = list(set(tuplesToInsert))
    
#    f = open('C:\Users\Simonial\Desktop\DRUGINFO.txt', 'w')
#    for i in tuplesToInsert:
#        f.write(str(i) + '\n')
#    f.close()
    
    # Insert the tuples into the database.
    values = '(' + ('%s,' * len(tuplesToInsert[0]))
    values = values[:-1] + ')'
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    mysql.tableINSERT(cursor, tableDrugs, values, tuplesToInsert)
    mysql.closeConnection(conn, cursor)
    
#    for i in targetDrugBinding.keys()[:10]:
#        print i, targetDrugBinding[i]