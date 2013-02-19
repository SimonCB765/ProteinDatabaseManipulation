'''
Created on 11 Oct 2011

@author: Simon Bull
'''

import re

import utilities.list2file
import utilities.XMLparser

def main(approvedTargetInputFile, XMLInputFile, DBTargetExternalLinks, targetOutputFile, drugOutputFile):
    """
    Takes three files containging the approved targets, drug information and external database cross-references for all targets. Returns two files.
        targetOutputFile contains the mapping of UniProt accessions to drugs that are approved and target it.
        drugOutputFile contains all the desired information from DrugBank about the individual drugs.
    targetOutputFile - A tab separated (tsv) file, with two elements on each line.
        The first element is a UniProt accession of an approved drug target.
        The second element is a semi-colon separated list of all the drugs that are approved and target the protein in the first element.
    drugOutputFile - A tab separated (tsv) file, with five elements on each line.
        The first element is the DrugBank ID of the drug.
        The second element is the name of the drug as recorded by DrugBank.
        The third element is a semi-colon separated list of all the DrugBank drug groups that the drug is a member of.
        The fourth element is the CAS number of the drug as recorded by DrugBank.
        The fifth element is a semi-colon separated list of PubChem CIDs that DrugBank has linked to the drugs.
    """

    # Extract drug information from the XML file of the DrugTarget database.
    allDrugsFromXML = []  # Records tuples containing the information of interest about every drug in DrugBank.
    approvedDrugIDs = set([])  # Records all of the approved drugs found in DrugBank.
    parser = utilities.XMLparser.XMLParser(XMLInputFile)
    parser.parse()
    allInfo = parser.retrieve(['drugs.drug', 'drugs.drug.drugbank-id', 'drugs.drug.name', 'drugs.drug.cas-number',
                               'drugs.drug.groups.group', 'drugs.drug.external-identifiers'])

    for i in allInfo:
        # Go through the record of each drug, recording the information of interest.
        drugType = allInfo[i]['drugs.drug'][0].parameters['type']
        if drugType == 'small molecule':
            # If the drug is a small molecule, then we are interested in it.
            drugID = allInfo[i]['drugs.drug.drugbank-id'][0].data  # Record the DrugBank ID of the drug.
            groups = ';'.join([j.data for j in allInfo[i]['drugs.drug.groups.group']])  # Determine the groups that the drug is member of.
            if 'approved' in groups:
                # If the drug is an approved drug, then add it the set of approved drugs.
                approvedDrugIDs.add(drugID)
            drugName = allInfo[i]['drugs.drug.name'][0].data  # Record the name of the drug.
            CASNumber = allInfo[i]['drugs.drug.cas-number'][0].data  # Record the CAS number of the drug.
            PubChemCIDs = []  # Used to record all the PubChem CIDs that are recorded for the drug.
            if allInfo[i]['drugs.drug.external-identifiers'] != []:
                for j in allInfo[i]['drugs.drug.external-identifiers'][0].children:
                    if j.children[0].data == 'PubChem Compound':
                        # If the data containined in the <drugs><drug><external-identifiers><external-identifier><resource> tag is 'PubChem Compound',
                        # then a PubChem CID is containined in the corresponding <drugs><drug><external-identifiers><external-identifier><identifier> tag.
                        PubChemCIDs.append(j.children[1].data)
            PubChemCIDs = ';'.join(PubChemCIDs)
            allDrugsFromXML.append(drugID + '\t' + drugName + '\t' + groups + '\t' + CASNumber + '\t' + PubChemCIDs)

    # Read the CSV file of external target links, and determine the UniProt accession number of the DrugBank targets.
    DBToUP = {}  # Record a mapping of DrugBank target ID to UniProt accession.
    readLinks = open(DBTargetExternalLinks, 'r')
    for line in readLinks:
        result = re.findall(',".*?",', line)  # Some of the target names contain commas within the name (and are therefore enclosed in ""). Find them.
        for i in result:
            # As the commas in the target names messes up the CSV format, and I'm not interested in the target names here, replace all target names that have a comma in them with a blank.
            # this allows proper splitting of the CSV file into its comma delimited chunks.
            line = line.replace(i, ',,')
        chunks = line.split(',')
        DBToUP[chunks[0]] = chunks[5]
    readLinks.close()

    # Extract target and sequence information from the list of approved drug targets.
    approvedTargetsFromTarget = []  # Record all the approved drugs that target an approved drug target.

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
                # If no UniProt accession was recorded for the target, then ignore it as no cross-referencing with the UniProt protein properties can be performed.
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