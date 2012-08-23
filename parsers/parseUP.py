'''
Created on 12 Oct 2011

@author: Simon Bull
'''

import os
import shutil
import re

import utilities.XMLparser
import utilities.list2file

def main(XMLInputFile, UPGPCRs, UPKinases, UPProteases, UPHumanAccessions, UPHumanAccessionMap, UPHumanNames,
         UPDrugIDs, UPProteinInfo, UPExternalLinks, UPPPIData, folderUP):
    """
    """
    
    # Extract the UP accessions from the GPCR file.
    readFile = open(UPGPCRs, 'r')
    gpcrUPAccs = []
    for line in readFile:
        foundHuman = re.search('(?<=_)HUMAN', line)
        if foundHuman:
            foundUPID = re.search('(?<=\()[a-zA-Z][a-zA-Z0-9]{5}(?=\))', line)
            foundUPID = foundUPID.group(0)
            gpcrUPAccs.append(foundUPID)
    readFile.close()
    
    # Extract the UP accessions from the kinase file.
    readFile = open(UPKinases, 'r')
    kinaseUPAccs = []
    for line in readFile:
        foundHuman = re.search('(?<=_)HUMAN', line)
        if foundHuman:
            foundUPID = re.search('(?<=\()[a-zA-Z][a-zA-Z0-9]{5}(?=\))', line)
            foundUPID = foundUPID.group(0)
            kinaseUPAccs.append(foundUPID)
    readFile.close()
    
    # Extract the UP accessions from the protease file.
    readFile = open(UPProteases, 'r')
    proteaseUPAccs = []
    for line in readFile:
        foundHuman = re.search('_HUMAN \([a-zA-Z][a-zA-Z0-9]{5}\),', line)
        if foundHuman:
            line = foundHuman.group(0)
            foundUPID = re.search('(?<=\()[a-zA-Z][a-zA-Z0-9]{5}(?=\))', line)
            foundUPID = foundUPID.group(0)
            proteaseUPAccs.append(foundUPID)
    readFile.close()
    
    modeOfActionKeywords = {'Kinase' : set(['KW-0418', 'KW-0723', 'KW-0829']),
                            'GPCR' : set(['KW-0297']),
                            'IonChannel' : set(['KW-1071', 'KW-0851', 'KW-0107', 'KW-0869', 'KW-0407', 'KW-0631',
                                            'KW-0894']),
                            'Protease' : set(['KW-0031', 'KW-0064', 'KW-0121', 'KW-0224', 'KW-0482', 'KW-0645',
                                          'KW-0720', 'KW-0788', 'KW-0888'])
                            }
    
    UPRepresentativeAccessions = []
    UPProteinNames = []
    UPAccessionMapping = []
    UPDrugMapping = []
    proteinTuples = []
    proteinExternalLinks = []
    
    splitXMLLocations = split_uniprot(XMLInputFile, 200, folderUP)
    
    PPIData = []
    for XMLFile in splitXMLLocations:
        parser = utilities.XMLparser.XMLParser(XMLFile)
        parser.parse()
        
        tags = parser.retrieve(['uniprot.entry.accession', 'uniprot.entry.name',
                                'uniprot.entry.comment',
                                'uniprot.entry.keyword', 'uniprot.entry.dbReference', 'uniprot.entry.feature',
                                'uniprot.entry.sequence'])
        for i in tags.keys():
            proteinRecord = ['acc', 'name', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0, 0, 0.0, 0.0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                             'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'N', 'seq']
            databaseTypes = {'GO' : [], 'GeneID' : [], 'UniGene' : [], 'DrugBank' : [], 'EC' : [], 'HGNC' : []}
            
            accessions = [j.data for j in tags[i]['uniprot.entry.accession']]
            representativeAcc = accessions[0]
            UPRepresentativeAccessions.append(representativeAcc)
            proteinRecord[0] = representativeAcc
            for acc in accessions:
                UPAccessionMapping.append(acc + '\t' + representativeAcc)
            
            proteinName = tags[i]['uniprot.entry.name'][0].data
            proteinRecord[1] = proteinName
            UPProteinNames.append(proteinName)
            
            comments = [j for j in tags[i]['uniprot.entry.comment']]
            isoforms = []
            subcellularLocation = []
            for j in comments:
                if j.parameters['type'] == 'subcellular location':
                    if j.children[0].tagName == 'molecule':
                        isoformLocation = j.children[0].data[8:]  # Isoform info is recorded as 'Isoform XXX', so remove the 'Isoform ' from the front.
                    else:
                        isoformLocation = 'MatureProtein'
                    for k in j.children:
                        if k.tagName == 'subcellularLocation' and k.children[0].tagName == 'location':
                            subcellularLocation.append(isoformLocation + ',' + k.children[0].data)
                elif j.parameters['type'] == 'alternative products':
                    for k in j.children:
                        if k.tagName == 'isoform':
                            isoforms.append(k.children[0].data + ',' + k.children[1].data)
                elif j.parameters['type'] == 'interaction':
                    proteinOne = j.children[0].parameters['intactId']
                    proteinTwo = j.children[1].parameters['intactId']
                    if proteinOne != proteinTwo:
                        # Not homotypic (protein not interacting with itself).
                        proteinTwo = j.children[1].children[0].data
                        if proteinTwo.count('-') != 0:
                            # If the protein is an isoform.
                            isoChunks = proteinTwo.split('-')
                            proteinTwo = isoChunks[0]
                            isoformID = isoChunks[1]
                        else:
                            isoformID = 'No Isoform'
                    else:
                        proteinTwo = representativeAcc
                        isoformID = 'No Isoform'
                    organismsDiffer = j.children[2].data
                    numExperiments = j.children[3].data
                    PPIData.append(tuple([representativeAcc, proteinTwo, isoformID, organismsDiffer, numExperiments]))
            subcellularLocation = ';'.join(subcellularLocation)
            if subcellularLocation == '':
                proteinRecord[43] = 'NA'
            else:
                proteinRecord[43] = subcellularLocation
            isoforms = ';'.join(isoforms)
            if isoforms == '':
                proteinRecord[53] = 'NA'
            else:
                proteinRecord[53] = isoforms
            
            modes = [j.parameters['id'] for j in tags[i]['uniprot.entry.keyword']]
            if modeOfActionKeywords['GPCR'].intersection(modes) != set([]):
                modeOfAction = 'G-protein coupled receptor'
            elif modeOfActionKeywords['IonChannel'].intersection(modes) != set([]):
                modeOfAction = 'Ion Channel'
            elif modeOfActionKeywords['Kinase'].intersection(modes) != set([]):
                modeOfAction = 'Kinase'
            elif modeOfActionKeywords['Protease'].intersection(modes) != set([]):
                modeOfAction = 'Protease'
            else:
                modeOfAction = 'NA'
            # Now that the modes extracted from the XML file have been used, annotate with the modes extracted from
            # the gpcr, kinase and protease files.
            if representativeAcc in gpcrUPAccs:
                modeOfAction = 'G-protein coupled receptor'
            elif representativeAcc in kinaseUPAccs:
                modeOfAction = 'Kinase'
            elif representativeAcc in proteaseUPAccs:
                modeOfAction = 'Protease'
            proteinRecord[36] = modeOfAction
            
            databaseTags = [j for j in tags[i]['uniprot.entry.dbReference'] if j.parameters['type'] in databaseTypes.keys()]
            for j in databaseTags:
                type = j.parameters['type']
                if type == 'GeneID':
                    databaseTypes['GeneID'].append(j.parameters['id'])
                elif type == 'UniGene':
                    databaseTypes['UniGene'].append(j.parameters['id'])
                elif type == 'DrugBank':
                    databaseTypes['DrugBank'].append(j.parameters['id'])
                elif type == 'GO':
                    databaseTypes['GO'].append(j.parameters['id'])
                elif type == 'HGNC':
                    databaseTypes['HGNC'].append(j.parameters['id'].split(':')[1])  # Have to split as we only want the number bit of the HGNC ID.
#                elif type == 'Ensembl':
#                    gene = ''
#                    for k in j.children:
#                        if k.parameters['type'] == 'gene ID':
#                            gene = k.parameters['value']
#                    if gene not in databaseTypes['Ensembl']:
#                        # It is possible that UniProt has multiple Ensembl transcript and protein entries for the same
#                        # Ensembl gene entry. Only record the gene when it has not already been recorded.
#                        databaseTypes['Ensembl'].append(gene)
                elif type == 'EC':
                    databaseTypes['EC'].append(j.parameters['id'])
            externalLinks = (representativeAcc + ',' + ';'.join(databaseTypes['GeneID']) + ',' +
                             ';'.join(databaseTypes['UniGene']) + ',' + ';'.join(databaseTypes['GO']) + ',' + ';'.join(databaseTypes['HGNC']))
#            externalLinks = (representativeAcc + ',' + ';'.join(databaseTypes['GeneID']) + ',' +
#                             ';'.join(databaseTypes['UniGene']) + ',' + ';'.join(databaseTypes['GO']) + ',' +
#                             ';'.join(databaseTypes['Ensembl']))
            proteinExternalLinks.append(externalLinks)
            if databaseTypes['DrugBank'] != []:
                UPDrugMapping.append(representativeAcc + '\t' + ';'.join(databaseTypes['DrugBank']))
            if databaseTypes['EC'] == []:
                proteinRecord[37] = 'NA'
            else:
                proteinRecord[37] = ';'.join(databaseTypes['EC'])
            
            glycosylationO = []
            glycosylationN = []
            phosphSer = []
            phosphThr = []
            phosphTyr = []
            topologicalDomain = []
            signalPeptide = []
            transMem = []
            turns = []
            helices = []
            strands = []
            for j in tags[i]['uniprot.entry.feature']:
                params = j.parameters
                if params['type'] == 'glycosylation site':
                    locationTag = j.children[0]
                    position = locationTag.children[0].parameters['position']
                    glycType = params['description']
                    if glycType[:8] == 'N-linked':
                        glycosylationN.append(position)
                    elif glycType[:8] == 'O-linked':
                        glycosylationO.append(position)
                elif params['type'] == 'modified residue':
                    locationTag = j.children[0]
                    position = locationTag.children[0].parameters['position']
                    phosphType = params['description']
                    if phosphType[:13] == 'Phosphoserine':
                        phosphSer.append(position)
                    elif phosphType[:16] == 'Phosphothreonine':
                        phosphThr.append(position)
                    elif phosphType[:15] == 'Phosphotyrosine':
                        phosphTyr.append(position)
                elif params['type'] == 'topological domain':
                    domain = params['description']
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    topologicalDomain.append(domain + ',' + beginLoc + ',' + endLoc)
                elif params['type'] == 'signal peptide':
                    if params.has_key('evidence'):
                        # The signal peptide location has been experimentally proven
                        evidence = 'Y'
                    else:
                        evidence = 'N'
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    signalPeptide.append(beginLoc + ',' + endLoc + ',' + evidence)
                elif (params['type'] == 'transmembrane region' and params.has_key('description') and
                      params['description'][:7] == 'Helical'):
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    transMem.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'turn':
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    turns.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'helix':
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    helices.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'strand':
                    locationTag = j.children[0]
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        endLoc = ''
                    strands.append(beginLoc + ',' + endLoc)
             
            proteinRecord[38] = ';'.join(glycosylationO) if glycosylationO != [] else 'NA'
            proteinRecord[39] = ';'.join(glycosylationN) if glycosylationN != [] else 'NA'
            proteinRecord[40] = ';'.join(phosphSer) if phosphSer != [] else 'NA'
            proteinRecord[41] = ';'.join(phosphThr) if phosphThr != [] else 'NA'
            proteinRecord[42] = ';'.join(phosphTyr) if phosphTyr != [] else 'NA'
            proteinRecord[44] = ';'.join(topologicalDomain) if topologicalDomain != [] else 'NA'
            proteinRecord[46] = ';'.join(signalPeptide) if signalPeptide != [] else 'NA'
            proteinRecord[47] = ';'.join(transMem) if transMem != [] else 'NA'
            proteinRecord[48] = ';'.join(turns) if turns != [] else 'NA'
            proteinRecord[49] = ';'.join(helices) if helices != [] else 'NA'
            proteinRecord[50] = ';'.join(strands) if strands != [] else 'NA'
#            tuple = (tuple + glycosylationO + '\t' + glycosylationN + '\t' + phosphSer + '\t' + phosphThr + '\t' +
#                     phosphTyr + '\t' + signalPeptide + '\t' + transMem + '\t' + turns + '\t' + helices + '\t' +
#                     strands + '\t')
    
            sequence = (tags[i]['uniprot.entry.sequence'][0].data).replace(' ', '')
            proteinRecord[-1] = sequence
            
            proteinTuples.append(tuple(proteinRecord))
    
    utilities.list2file.main(UPRepresentativeAccessions, UPHumanAccessions)
    utilities.list2file.main(UPAccessionMapping, UPHumanAccessionMap)
    utilities.list2file.main(UPProteinNames, UPHumanNames)
    utilities.list2file.main(UPDrugMapping, UPDrugIDs)
    utilities.list2file.main(proteinTuples, UPProteinInfo)
    utilities.list2file.main(proteinExternalLinks, UPExternalLinks)
    utilities.list2file.main(PPIData, UPPPIData)

def split_uniprot(XMLInputFile, entriesPerSplit, folderUP):
    
    splitXMLLocations = []
    
    splitDir = folderUP + '/TempXML'
    if os.path.exists(splitDir):
        shutil.rmtree(splitDir)
    os.mkdir(splitDir)
    
    entriesInCurrentSplit = 0
    newFile = True
    currentFileNumber = 0
    readIn = open(XMLInputFile, 'r')
    for line in readIn:
        if newFile:
            if line == '</uniprot>\n':
                break
            currentOutputFile = splitDir + '/File' + str(currentFileNumber) + '.xml'
            writeOut = open(currentOutputFile, 'w')
            newFile = False
            if currentFileNumber != 0:
                writeOut.write('<uniprot>\n')
        
        writeOut.write(line)
        
        if line == '</entry>\n':
            entriesInCurrentSplit += 1
        if entriesInCurrentSplit == entriesPerSplit:
            writeOut.write('</uniprot>\n')
            writeOut.close()
            currentFileNumber += 1
            entriesInCurrentSplit = 0
            splitXMLLocations.append(currentOutputFile)
            newFile = True
    try:
        writeOut.close()
        splitXMLLocations.append(currentOutputFile)
    except:
        pass
    
    return splitXMLLocations