'''
Created on 16 Nov 2011

@author: Simon Bull
'''

def main(bindingSDF, binding2PubChem, bindingParsed):
    """
    Takes two files containing information about interactions between compounds and targets and information about BindingDB compound ID to PubChem CID mappings.
    Returns a file of binding constants between compounds and proteins.
    bindingParsed - A tab separated file (tsv), with four elements on each line.
        The first element is the UniProt accessions targeted by a compound.
        The second element is the CID of the compound that targets the UniProt accessions in the first element.
        The third element is the Ki for the interaction between the compound in element two and the proteins in element one.
        The fourth element is the Kd for the interaction between the compound in element two and the proteins in element one.
    """

    # Determine the mapping of BindingDB compound IDs to PubChem CIDs.
    pubChemDict = {}
    readPubChem = open(binding2PubChem, 'r')
    for line in readPubChem:
        chunks = line.split()
        pubChemDict[chunks[0]] = chunks[1]
    readPubChem.close()

    moleculeDict = {}
    pubChemIDsUsed = set([])
    bindingID = 0
    targetNumber = ''
    targetUP = ''
    affinityMeasure = ''
    affinityUnits = ''
    humanTarget = False  # Records whether the target is human.
    recordBindingDBID = False  # Records whether you are expecting to see the compound ID on the next line.
    recordOrganism = False  # Records whether you are expecting to find the target organism on the next line.
    recordUP = False  # Records whether you are expecting to see the UniProt accession of the target on the next line.
    recordBindingAffinity = False  # Records whether you are expecting measurement (Kd, Ki, etc.) information on the next line.

    readFrom = open(bindingSDF, 'r')
    writeTo = open(bindingParsed, 'w')
    for line in readFrom:
        if  line == '$$$$\n':
            # If you've reached the end of a molecule entry, record the information.
            uniqueRecords = set([])
            if pubChemDict.has_key(bindingID) and moleculeDict != {}:
                for i in moleculeDict.keys():
                    Ki = moleculeDict[i]['Ki'][0]
                    Kd = moleculeDict[i]['Kd'][0]
                    uniqueRecords.add(i + '\t' + pubChemDict[bindingID] + '\t' + Ki + '\t' + Kd)
            for i in uniqueRecords:
                writeTo.write(str(i) +'\n')

        if recordBindingDBID:
            # If you are expecting to be recording the BindingDB compound ID on this line.
            moleculeDict = {}  # Initialise a new empty dictionary to hold the binding data about the compound.
            bindingID = line.strip()  # Record the compound ID.
            recordBindingDBID = False  # Indicate that you are not expecting to find the compound ID on the next line.
        elif recordOrganism:
            # If you are expecting to find the organism of a target on this line.
            if line.strip() == 'Homo sapiens':
                # If the organism is human, then indicate this.
                humanTarget = True
            else:
                humanTarget = False
            recordOrganism = False  # Indicate that you are not expecting to find the organism type on the next line.
        elif humanTarget:
            # Only record measurement and UniProt accessions if the target is human.
            if recordUP:
                targetUP = ';'.join(line.split())  # There may be more than one UniProt accession for the target, so turn this into a semi-colon separated list of accessons.
                if not moleculeDict.has_key(targetUP):
                    # If the target UniProt accession has not been seen before, then initiliase it with no known tightest binding constant information.
                    moleculeDict[targetUP] = {'Ki' : [], 'Kd' : []}
                recordUP = False  # Indicate that you are not expecting to find a UniProt accession on the next line.
            elif recordBindingAffinity:
                value = line.strip()
                if affinityMeasure == 'Ki':
                    moleculeDict[targetUP][affinityMeasure].append(value)
                elif affinityMeasure == 'Kd':
                    moleculeDict[targetUP][affinityMeasure].append(value)
                recordBindingAffinity = False  # Indicate that you are not expecting to find measurement information onthe next line.
        else:
            recordUP = False  # Indicate that you are not expecting to find a UniProt accession on the next line.
            recordBindingAffinity = False  # Indicate that you are not expecting to find measurement information onthe next line.

        if line == '> <BindingDB monomerid>\n':
            # Indicate that the next line will have the ID of the compound on it.
            recordBindingDBID = True
        elif line[:25] == '> <TARGET Source Organism':
            # Indicate that the next line will have the target organism.
            recordOrganism = True
        elif line[:29] == '> <UniProtKB Accession Number':
            # Indicate that the next line will have the UniProt accession of the target on it.
            recordUP = True
        elif line[:15] == '> <Enzymologic:':
            chunks = line.split()
            affinityMeasure = chunks[2]  #Record the type of measurement (Ki, Kd, etc.) is being made.
            recordBindingAffinity = True  # Indicate that the next line is expected to have measurement information on it.

    writeTo.close()
    readFrom.close()