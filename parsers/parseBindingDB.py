'''
Created on 16 Nov 2011

@author: Simon Bull
'''

#import gzip
#
#import utilities.PubChemSDFdownload
#import utilities.file2list
#import utilities.list2file

def main(bindingSDF, binding2PubChem, bindingParsed):#, bindingPubChemCIDs, bindingPubChemSDF, bindingInChis):
    
    # Parse the BindingDB SDF file.
    print '\tParsing BindingDB SDF file.'
    parse_binding_sdf(bindingSDF, binding2PubChem, bindingParsed)#, bindingPubChemCIDs)
    
#    # Get a gzipped sdf file containing the PubChem compounds linked to the BindingDB compounds (http://pubchem.ncbi.nlm.nih.gov/pug/pughelp.html#31).
#    print '\tDownloading structures from PubChem.'
#    listOfCIDs = utilities.file2list.main(bindingPubChemCIDs)
#    listOfCIDs = [eval(i) for i in listOfCIDs]
#    
#    blockSize = 10000  # There is a limit of 250,000 structure per query.
#    lowerBound = 0
#    upperBound = blockSize
#    inchiList = set([])
#    while lowerBound < len(listOfCIDs):
#        utilities.PubChemSDFdownload.main(listOfCIDs[lowerBound:upperBound], bindingPubChemSDF)
#        readFrom = open(bindingPubChemSDF, 'rb')
#        zipped = gzip.GzipFile(fileobj=readFrom)
#        recordCompoundID = False
#        recordInChi = False
#        currentRecord = ''
#        for line in zipped:
#            if recordCompoundID:
#                currentRecord += line.strip() + '\t'
#                recordCompoundID = False
#            elif recordInChi:
#                currentRecord += line.strip()
#                inchiList.add(currentRecord)
#                currentRecord = ''
#                recordInChi = False
#    
#            if line == '> <PUBCHEM_COMPOUND_CID>\n':
#                recordCompoundID = True
#            elif line == '> <PUBCHEM_IUPAC_INCHI>\n':
#                recordInChi = True
#        zipped.close()
#        readFrom.close()
#        lowerBound += blockSize
#        upperBound += blockSize
#    
#    utilities.list2file.main(list(inchiList), bindingInChis)

def parse_binding_sdf(bindingSDF, binding2PubChem, bindingParsed):#, bindingPubChemCIDs):
    
    pubChemDict = {}
    readPubChem = open(binding2PubChem, 'r')
    for line in readPubChem:
        chunks = line.split()
        pubChemDict[chunks[0]] = chunks[1]
    readPubChem.close()
    
    readFrom = open(bindingSDF, 'r')
    writeTo = open(bindingParsed, 'w')
    
    moleculeDict = {}
    pubChemIDsUsed = set([])
    
    bindingID = 0
    targetNumber = ''
    targetUP = ''
    affinityMeasure = ''
    affinityUnits = ''
    
    humanTarget = False
    recordBindingDBID = False
    recordOrganism = False
    recordUP = False
    recordBindingAffinity = False
    
    for line in readFrom:
        if  line == '$$$$\n':
            # If you've reached the end of a molecule entry, record the information.
            uniqueRecords = set([])
            if pubChemDict.has_key(bindingID) and moleculeDict != {}:
                for i in moleculeDict.keys():
                    Ki = moleculeDict[i]['Ki'][0]
##                    KiUnits = moleculeDict[i]['Ki'][1]
##                    IC50 = moleculeDict[bindingID]['Targets'][i]['IC50'][0]
##                    IC50Units = moleculeDict[bindingID]['Targets'][i]['IC50'][1]
                    Kd = moleculeDict[i]['Kd'][0]
##                    KdUnits = moleculeDict[i]['Kd'][1]
##                    EC50 = moleculeDict[bindingID]['Targets'][i]['EC50/IC50'][0]
##                    EC50Units = moleculeDict[bindingID]['Targets'][i]['EC50/IC50'][1]
                    uniqueRecords.add(i + '\t' + pubChemDict[bindingID] + '\t' + Ki + '\t' + Kd)
            for i in uniqueRecords:
                writeTo.write(str(i) +'\n')
    
        if recordBindingDBID:        
            moleculeDict = {}
            bindingID = line.strip()
            recordBindingDBID = False
        elif recordOrganism:
            if line.strip() == 'Homo sapiens':
                humanTarget = True
##                moleculeDict[bindingID]['Targets'][targetNumber] = {}
            else:
                humanTarget = False
            recordOrganism = False
        elif humanTarget:
            if recordUP:
                targetUP = ';'.join(line.split())
                if not moleculeDict.has_key(targetUP):
                    moleculeDict[targetUP] = {'Ki' : [], 'Kd' : []}
                recordUP = False
            elif recordBindingAffinity:
                value = line.strip()
                if affinityMeasure == 'Ki':
                    moleculeDict[targetUP][affinityMeasure].append(value)
##                elif affinityMeasure == 'IC50':
##                    moleculeDict[bindingID]['Targets'][targetNumber][affinityMeasure] = [value, affinityUnits]
                elif affinityMeasure == 'Kd':
                    moleculeDict[targetUP][affinityMeasure].append(value)
##                elif affinityMeasure == 'EC50/IC50':
##                    moleculeDict[bindingID]['Targets'][targetNumber][affinityMeasure] = [value, affinityUnits]
                recordBindingAffinity = False
        else:
            recordUP = False
            recordBindingAffinity = False
    
        if line == '> <BindingDB monomerid>\n':
            recordBindingDBID = True
        elif line[:21] == '> <TARGET Biomolecule':
            humanTarget = True
            chunks = (line.strip()).split()
##            targetNumber = chunks[-1][:-1]
        elif line[:25] == '> <TARGET Source Organism':
            recordOrganism = True
        elif line[:29] == '> <UniProtKB Accession Number':
            recordUP = True
        elif line[:15] == '> <Enzymologic:':
            chunks = line.split()
            affinityMeasure = chunks[2]
##            affinityUnits = chunks[3]
            recordBindingAffinity = True
    
    writeTo.close()
    readFrom.close()