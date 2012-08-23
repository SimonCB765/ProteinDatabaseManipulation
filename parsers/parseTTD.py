'''
Created on 11 Oct 2011

@author: Simon Bull
'''

import re

import utilities.list2file

def main(TTDTargets, TTDUPAccessions, TTDDrugXref, TTDTarget2Drug, TTDCIDs):
    """Extract only 'successful' (i.e. approved) targets from the TTD database.

    The database should be downloaded in the raw text format.

    TTDTargets @type - string
    TTDTargets @use  - the file location of the TTD database
    TTDUPAccessions @type - string
    TTDUPAccessions @use  - the file location for the UniProt Accessions of the approved TTD targets to be output
    """
    
    foundProtIDs = []
    target2Drug = {}
    target2Uniprot = {}

    inputFile = open(TTDTargets,'r')
    for line in inputFile:
        if line[:4] == 'TTDS':
            # Found an approved target.
            
            # It's possible for one TTD target to be associated with multiple UniProt IDs.
            # For this reason all UniProt IDs associated with a TTD target are recorded.
            chunks = (line.strip()).split('\t')
            if chunks[1] == 'UniProt ID':
                # If this is true the line contains a UniProt ID.
                foundProtIDs.append(chunks[2])
                if target2Uniprot.has_key(chunks[0]):
                    target2Uniprot[chunks[0]].append(chunks[2])
                else:
                    target2Uniprot[chunks[0]] = []
                    target2Uniprot[chunks[0]].append(chunks[2])
            elif chunks[1] == 'Drug(s)' and chunks[5] == 'Approved':
                # If the line contains an approved drug.
                if target2Drug.has_key(chunks[0]):
                    target2Drug[chunks[0]].append([chunks[2], chunks[3]])
                else:
                    target2Drug[chunks[0]] = []
                    target2Drug[chunks[0]].append([chunks[2], chunks[3]])
    inputFile.close()

    foundProtIDs.sort()

    utilities.list2file.main(foundProtIDs, TTDUPAccessions)
    
    # Extract drug info.
    pubchemCIDs = set([])
    drugInfo = {}
    inputFile = open(TTDDrugXref, 'r')
    for line in inputFile:
        chunks = (line.strip()).split('\t')
        if not re.search('[A-Z]{3}[0-9]{6}', chunks[0]):
            continue
        if chunks[2] == 'CAS Number':
            drugInfo[chunks[0]] = {'CAS' : '', 'CID' : ''}
            CAS = chunks[3].split()
            drugInfo[chunks[0]]['CAS'] = CAS[1]
        elif chunks[2] == 'PubChem CID':
            drugCID = chunks[3].split()
            pubchemCIDs.add(drugCID[1])
            if not drugInfo.has_key(chunks[0]):
                drugInfo[chunks[0]] = {'CAS' : '', 'CID' : ''}
                drugInfo[chunks[0]]['CID'] = drugCID[1]
            else:
                drugInfo[chunks[0]]['CID'] = drugCID[1]
    inputFile.close()
    
    target2DrugOutput = []
    for i in target2Drug.keys():
        target = i
        if not target2Uniprot.has_key(i):
            # Only bother with targets that have UniProt Accessions, as there is no way to cross reference targets
            # without them
            continue
        uniprot = ','.join(target2Uniprot[i])
        drugs = []
        for j in target2Drug[i]:
            drugName = j[0]
            drugID = j[1]
            if drugInfo.has_key(drugID):
                CAS = drugInfo[drugID]['CAS']
                if drugInfo[drugID]['CID'] == '':
                    CID = ''
                else:
                    CID = drugInfo[drugID]['CID']
            drugs.append(drugName + ',' + CAS + ',' + CID)
        drugs = ';'.join(drugs)
        target2DrugOutput.append('\t'.join([target, uniprot, drugs]))
    
    utilities.list2file.main(target2DrugOutput, TTDTarget2Drug)
    
    print '\tNumber of approved UniProt targets in the TTD: ', len(target2DrugOutput)