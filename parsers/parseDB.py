'''
Created on 11 Oct 2011

@author: Simon Bull
'''

import re

import utilities.list2file
import utilities.XMLparser

def main(approvedTargetInputFile, XMLInputFile, DBTargetExternalLinks, targetOutputFile, drugOutputFile):
    """Used to extract DrugBank approved drug IDs and DrugBank (target ID, target sequence) pairs from DrugBank downloads.

    The DrugBank drug IDs are extracted from the XML file of the entire DrugBank database.
    The target IDs are extracted from the FASTA format file of the DrugBank approved targets.
    
    @param approvedTargetInputFile: the pathname of the location of the file to extract target information from
    @type approvedTargetInputFile:  string
    @param XMLInputFile: the pathname of the location of the file to extract drug information from
    @type XMLInputFile:  string
    @param targetOutputFile: the pathname of the location to write the list of (target ID, target sequence) pairs to
    @type targetOutputFile:  string
    @param drugOutputFile: the pathname of the location to write the list of approved drug IDs to
    @type drugOutputFile:  string
    """
    
    # Extract drug information from the XML file of the DrugTarget database.
    allDrugsFromXML = []
    approvedDrugIDs = set([])
    parser = utilities.XMLparser.XMLParser(XMLInputFile)
    parser.parse()
    allInfo = parser.retrieve(['drugs.drug', 'drugs.drug.drugbank-id', 'drugs.drug.name', 'drugs.drug.cas-number',
                               'drugs.drug.groups.group', 'drugs.drug.external-identifiers'])
    for i in allInfo:
        drugType = allInfo[i]['drugs.drug'][0].parameters['type']
        if drugType == 'small molecule':
            drugID = allInfo[i]['drugs.drug.drugbank-id'][0].data
            groups = ';'.join([j.data for j in allInfo[i]['drugs.drug.groups.group']])
            if 'approved' in groups:
                approvedDrugIDs.add(drugID)
            drugName = allInfo[i]['drugs.drug.name'][0].data
            CASNumber = allInfo[i]['drugs.drug.cas-number'][0].data
            PubChemCIDs = []
            if allInfo[i]['drugs.drug.external-identifiers'] != []:
                for j in allInfo[i]['drugs.drug.external-identifiers'][0].children:
                    if j.children[0].data == 'PubChem Compound':
                        PubChemCIDs.append(j.children[1].data)
            PubChemCIDs = ';'.join(PubChemCIDs)
            allDrugsFromXML.append(drugID + '\t' + drugName + '\t' + groups + '\t' + CASNumber + '\t' + PubChemCIDs)   
    
    # Read the CSV file of external target links, and determine the UniProt accession number of the DrugBank targets.
    DBToUP = {}
    readLinks = open(DBTargetExternalLinks, 'r')
    for line in readLinks:
        result = re.findall(',".*?",', line)
        for i in result:
            line = line.replace(i, ',,')
        chunks = line.split(',')
        DBToUP[chunks[0]] = chunks[5]
    readLinks.close()
    
    # Extract target and sequence information from the list of approved drug targets.
    approvedTargetsFromTarget = []
    
    readTargets = open(approvedTargetInputFile, 'r')
    
    for line in readTargets:
        if line[0] == '>':
            # Found the start of a protein entry.
            # Split the line into two parts with the part containing the target ID at chunks[0].
            chunks = line.split(' ', 1)
            targetID = chunks[0]
            # Split the target ID into two with the numeric portion of the ID in targetID[1].
            targetID = targetID.split('|')
            targetID = targetID[1]
            if not DBToUP.has_key(targetID):
                continue
            targetID = DBToUP[targetID]
            # Get the DrugBank drug IDs of the drugs that target the target.
            drugs = re.findall('DB[0-9]{5}', line)
            drugs = [i for i in drugs if i in approvedDrugIDs]
            if drugs == []:
                # If there are no approved drugs that target the protein (caused by all approved drugs targetting it
                # being biotech drugs)
                continue
            drugs = ';'.join(drugs)
            approvedTargetsFromTarget.append(targetID + '\t' + drugs)
    readTargets.close()
    
    allDrugs = list(set(allDrugsFromXML))
    print '\tNumber of approved drugs in DrugBank:', len(allDrugs)
    # Provided you have found some drugs write the information out
    utilities.list2file.main(allDrugs, drugOutputFile)

    print '\tNumber of approved targets in DrugBank:', len(approvedTargetsFromTarget)
    # Provided you have found some targets write the information out
    utilities.list2file.main(approvedTargetsFromTarget, targetOutputFile)