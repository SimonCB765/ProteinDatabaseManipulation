'''
Created on 11 Oct 2011

@author: Simon Bull
'''

import re

import utilities.list2file

def main(TTDTargets, TTDUPAccessions, TTDDrugXref, TTDTarget2Drug):
    """
    Takes the TTD target database data and the TTD drug cross-reference data, and returns ont file of the UniProt accessions of the approved target proteins and one containing a mapping of approved targets to approved drugs.
    TDUPAccessions - A file of UniProt accessions, one on each line.
        Each UP accession in the file is the target of an approved drug (as recorded by the TTD).
    TTDTarget2Drug - A tab separated (tsv) file, with three elements on each line.
        The first element is the TTD ID for the protein.
        The second element is a comma separated list of UniProt accessions that the first element is linked to.
        The third element is a semi-colon separated list of approved drugs that target the protein. For each drug a comma separated list of three elements is recorded.
            The first element is the name of the drug (as recorded by the TTD).
            The second element is the CAS number of the drug.
            The third element is the PubChem CID for the drug.
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