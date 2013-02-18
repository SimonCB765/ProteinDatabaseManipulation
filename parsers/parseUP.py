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
    Takes the xml file of UniProt human proteins; files of UniProt GPCRs, kinases and proteases
    UPHumanAccessions - A file with one representative human protein UniProt accession on each line.
        The number of accessions in this file is the same as the number of human proteins recorded in UniProt.
    UPHumanAccessionMap - A tab separated (tsv) file, with two elements on each line.
        The first element is a UniProt protein accession (representative or not (most often not)).
        The second elements is the representative UniProt accession that the accession in the first element maps to.
    UPHumanNames - A file with one human protein name on each line.
        The number of names in this file is the same as the number of human proteins recorded in UniProt.
    UPDrugIDs - A tab separated (tsv) file, with two elements on each line.
        The first element is a representative UniProt accession.
        The second element is a semi-colon separated list of DrugBank drug IDs, one ID for each DrugBank drug that is recorded as being linked to the accession in UniProt.
    UPProteinInfo - A file of 56-tuples, one on each line.
        See the README for a full description of the file.
    UPExternalLinks - A comma separated (csv) file, with five (six after parseEnsembl is used) elements on each line.
        The first element is a representative UniProt accession.
        The second element is a semi-colon separated list of the Entrez Gene IDs that are recorded as being linked to the accession in UniProt.
        The third element is a semi-colon separated list of the UniGene cluster IDs that are recorded as being linked to the accession in UniProt.
        The fourth element is a semi-colon separated list of the Gene Ontology term IDs that are recorded as being linked to the accession in UniProt.
        The fifth element is a semi-colon separated list of the HGNC IDs that are recorded as being linked to the accession in UniProt.
        The sixth element is a semi-colon separated list of '-' separated lists of Ensembl records that are recorded as being linked to the accession in UniProt. The format for the '-' separated lists is:
            The first element is the Ensembl Gene ID that is recorded as being linked to the UniProt accession (because the transcript in the second element is linked to it).
            The second element is the Ensembl Transcript ID that is recorded as being linked to the UniProt accession (because the protein in the third element is linked to it).
            The third element is the Ensembl Protein ID that is recorded as being linked to the UniProt accession.
            Example : ENSG00000143627-ENST00000342741-ENSP00000339933;ENSG00000143627-ENST00000271946-ENSP00000271946
    UPPPIData - A file of 5-tuples, one on each line.
        The first element is the UniProt accession for the first protein in the interaction.
        The second element is the UniProt accession for the second protein in the interaction.
        The third element is a number if the second protein is an isofom (the number being the portion after the '-' in a UniProt isoform accession (e.g. 3 in O00257-3)), or 'No Isoform' if the second protein is not an isoform.
        The fourth element is false if the organism the second protein comes from is the same as that of the first (i.e. both are human), or true if the second protein is non-human.
        The fifth element is the number of experiments that give evidence for the interaction between the first and second proteins.
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

    # Define the UniProt keyword identifiers that indicate that the protein belongs to one of the four categories of interest.
    modeOfActionKeywords = {'Kinase' : set(['KW-0418', 'KW-0723', 'KW-0829']),
                            'GPCR' : set(['KW-0297']),
                            'IonChannel' : set(['KW-1071', 'KW-0851', 'KW-0107', 'KW-0869', 'KW-0407', 'KW-0631',
                                            'KW-0894']),
                            'Protease' : set(['KW-0031', 'KW-0064', 'KW-0121', 'KW-0224', 'KW-0482', 'KW-0645',
                                          'KW-0720', 'KW-0788', 'KW-0888'])
                            }

    UPRepresentativeAccessions = []  # Records the represetnative UniProt accession for every protein.
    UPProteinNames = []  # Records the name of every protein.
    UPAccessionMapping = []  # Records the mapping from every accession to its represetnative accession.
    UPDrugMapping = []  # Records a mapping from every represetnative protein accession to the drugs (approved or not) that target it
    proteinTuples = []  # Records the information about every represetnative accession that is to be recorded.
    proteinExternalLinks = []  # Records the external database identifiers linked to each represetnative accession.
    PPIData = []  # Records the information about which other proteins each representative accession participates in a binary PPI with.

    splitXMLLocations = split_uniprot(XMLInputFile, 200, folderUP)  # Split the XML file into manageable smaller files with 200 proteins per file.

    for XMLFile in splitXMLLocations:
        # For each of the smaller XML files.
        parser = utilities.XMLparser.XMLParser(XMLFile)
        parser.parse()

        # Extract the information of interest from the parsed XML files.
        tags = parser.retrieve(['uniprot.entry.accession', 'uniprot.entry.name',
                                'uniprot.entry.comment',
                                'uniprot.entry.keyword', 'uniprot.entry.dbReference', 'uniprot.entry.feature',
                                'uniprot.entry.sequence'])

        for i in tags.keys():
            # Define the default (empty) protein information tuple and database cross-reference link dictionary.
            proteinRecord = ['acc', 'name', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0, 0, 0.0, 0.0, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                             'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'N', 'seq']
            databaseTypes = {'GO' : [], 'GeneID' : [], 'UniGene' : [], 'DrugBank' : [], 'EC' : [], 'HGNC' : []}

            # Determine protein accession information.
            accessions = [j.data for j in tags[i]['uniprot.entry.accession']]  # Extract all the accessions for the current protein.
            representativeAcc = accessions[0]  # The represetnative accession is the first accession.
            UPRepresentativeAccessions.append(representativeAcc)
            proteinRecord[0] = representativeAcc
            for acc in accessions:
                # Add all non-representative to representative accession mappings for the current protein.
                UPAccessionMapping.append(acc + '\t' + representativeAcc)

            # Determine protein name information.
            proteinName = tags[i]['uniprot.entry.name'][0].data  # Extract the name of the protein.
            proteinRecord[1] = proteinName
            UPProteinNames.append(proteinName)

            # Determine the information of interest about the protein that resides in the comment section of the protein record.
            comments = [j for j in tags[i]['uniprot.entry.comment']]
            isoforms = []  # Used to record information about isoforms of the protein.
            subcellularLocation = []  # Used to record information about the subcellular locations that the protein is found in.
            for j in comments:
                if j.parameters['type'] == 'subcellular location':
                    # If the current comment is a record of subcellular location information.
                    if j.children[0].tagName == 'molecule':
                        # If the tag one level down from the <comment type="subcellular location" ...> tag is a <molecule> tag,
                        # then the subcellular location is for an isoform of the protein.
                        isoformLocation = j.children[0].data[8:]  # Isoform info is recorded as 'Isoform XXX', so remove the 'Isoform ' from the front.
                    else:
                        # The subcellular location information is for the protein itself (not an isoform).
                        isoformLocation = 'MatureProtein'
                    for k in j.children:
                        if k.tagName == 'subcellularLocation' and k.children[0].tagName == 'location':
                            # Record every subcellular location that the protein (or its isoform) occurs in.
                            # Ignore tags such as those indicating information about topogy or giving the evidence for the location.
                            subcellularLocation.append(isoformLocation + ',' + k.children[0].data)
                elif j.parameters['type'] == 'alternative products':
                    # The comment records alternative products of the protein (and potentially an isoform).
                    for k in j.children:
                        if k.tagName == 'isoform':
                            # The comment records information about an isoform.
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
                    organismsDiffer = j.children[2].data  # Determine whether the organisms differ (i.e. proteinTwo comes from a non-human source).
                    numExperiments = j.children[3].data  # Determine the number of experiments that give evidence for the interaction.
                    PPIData.append(tuple([representativeAcc, proteinTwo, isoformID, organismsDiffer, numExperiments]))
            subcellularLocation = ';'.join(subcellularLocation)
            if subcellularLocation == '':
                # If there was no subcellular information, then record it as NA.
                proteinRecord[43] = 'NA'
            else:
                proteinRecord[43] = subcellularLocation
            isoforms = ';'.join(isoforms)
            if isoforms == '':
                # If there are no isoforms for the protein, then record it as NA.
                proteinRecord[53] = 'NA'
            else:
                proteinRecord[53] = isoforms

            # Determine the information about the protein family that the protein is in (from the keyword section).
            modes = [j.parameters['id'] for j in tags[i]['uniprot.entry.keyword']]
            if modeOfActionKeywords['GPCR'].intersection(modes) != set([]):
                # If any of the GPCR keywords appear in the keywords for the protein, then the protein is a GPCR.
                modeOfAction = 'G-protein coupled receptor'
            elif modeOfActionKeywords['IonChannel'].intersection(modes) != set([]):
                # If any of the ion channel keywords appear in the keywords for the protein, then the protein is an ion channel.
                modeOfAction = 'Ion Channel'
            elif modeOfActionKeywords['Kinase'].intersection(modes) != set([]):
                # If any of the kinase keywords appear in the keywords for the protein, then the protein is a kinase.
                modeOfAction = 'Kinase'
            elif modeOfActionKeywords['Protease'].intersection(modes) != set([]):
                # If any of the protease keywords appear in the keywords for the protein, then the protein is a protease.
                modeOfAction = 'Protease'
            else:
                # If non of the protein family keywords appear in the keywords for the protein, then the protein is not in any of the four families of interest.
                modeOfAction = 'NA'
            # Now that the family information extracted from the XML file has been used, annotate with the family information extracted from the GPCR, kinase and protease files.
            if representativeAcc in gpcrUPAccs:
                # If the represetnative accession appears in the file of GPCRs, then the protein is a GPCR.
                modeOfAction = 'G-protein coupled receptor'
            elif representativeAcc in kinaseUPAccs:
                # If the represetnative accession appears in the file of kinases, then the protein is a kinase.
                modeOfAction = 'Kinase'
            elif representativeAcc in proteaseUPAccs:
                # If the represetnative accession appears in the file of proteases, then the protein is a protease.
                modeOfAction = 'Protease'
            proteinRecord[36] = modeOfAction

            # Extract the external database IDs that are linked to the current protein.
            databaseTags = [j for j in tags[i]['uniprot.entry.dbReference'] if j.parameters['type'] in databaseTypes.keys()]
            for j in databaseTags:
                type = j.parameters['type']
                if type == 'GeneID':
                    # The ID is an Entrez Gene ID.
                    databaseTypes['GeneID'].append(j.parameters['id'])
                elif type == 'UniGene':
                    # The ID is for a UniGene cluster.
                    databaseTypes['UniGene'].append(j.parameters['id'])
                elif type == 'DrugBank':
                    # The ID is for a DrugBank drug.
                    databaseTypes['DrugBank'].append(j.parameters['id'])
                elif type == 'GO':
                    # The ID is for a Gene Ontology term.
                    databaseTypes['GO'].append(j.parameters['id'])
                elif type == 'HGNC':
                    # The ID is for an HGNC gene ID.
                    databaseTypes['HGNC'].append(j.parameters['id'].split(':')[1])  # Have to split as we only want the number bit of the HGNC ID.
                elif type == 'EC':
                    # The ID is an EC number.
                    databaseTypes['EC'].append(j.parameters['id'])
            # Convert the external database IDs into the format expected for writing out.
            externalLinks = (representativeAcc + ',' + ';'.join(databaseTypes['GeneID']) + ',' +
                             ';'.join(databaseTypes['UniGene']) + ',' + ';'.join(databaseTypes['GO']) + ',' + ';'.join(databaseTypes['HGNC']))
            proteinExternalLinks.append(externalLinks)
            if databaseTypes['DrugBank'] != []:
                # If there are DrugBank drugs linked to the protein, then add this information to the record of accession DrugBank drug mappings.
                UPDrugMapping.append(representativeAcc + '\t' + ';'.join(databaseTypes['DrugBank']))
            if databaseTypes['EC'] == []:
                # If there was no EC number for the protein record this.
                proteinRecord[37] = 'NA'
            else:
                # Record the found EC number(s) for the protein.
                proteinRecord[37] = ';'.join(databaseTypes['EC'])

            glycosylationO = []  # Records the O-linked glycosylation sites for the current protein.
            glycosylationN = []  # Records the N-linked glycosylation sites for the current protein.
            phosphSer = []  # Records the serine phosphorylation sites for the current protein.
            phosphThr = []  # Records the threonine phosphorylation sites for the current protein.
            phosphTyr = []  # Records the tyrosine phosphorylation sites for the current protein.
            topologicalDomain = []  # Records the topological domains for the current protein.
            signalPeptide = []  # Records the signal peptide cleavage sites for the current protein.
            transMem = []  # Records the transmembrane regions for the current protein.
            turns = []  # Records the turns in the current protein.
            helices = []  # Records the alpha helices in the current protein.
            strands = []  # Records the beat strands in the current protein.
            for j in tags[i]['uniprot.entry.feature']:
                params = j.parameters
                if params['type'] == 'glycosylation site':
                    # The current feature is a glycosylation site.
                    locationTag = j.children[0]  # Extract the tag containing the location of the glycosylation site.
                    position = locationTag.children[0].parameters['position']  # Extract the site of the glycosylation.
                    glycType = params['description']
                    # Record whether the glycosylation is O-linked or N-linked.
                    if glycType[:8] == 'N-linked':
                        glycosylationN.append(position)
                    elif glycType[:8] == 'O-linked':
                        glycosylationO.append(position)
                elif params['type'] == 'modified residue':
                    # The current feature records modified residue information (and therefore contains the phosphorylation information we want).
                    locationTag = j.children[0]  # Get the tag containing the location of the modified residue.
                    position = locationTag.children[0].parameters['position']  # Get the site of the modified residue.
                    phosphType = params['description']
                    if phosphType[:13] == 'Phosphoserine':
                        # Record that the modified residue is a phosphorylated serine.
                        phosphSer.append(position)
                    elif phosphType[:16] == 'Phosphothreonine':
                        # Record that the modified residue is a phosphorylated threonine/
                        phosphThr.append(position)
                    elif phosphType[:15] == 'Phosphotyrosine':
                        # Record that the modified residue is a phosphorylated tyrosine.
                        phosphTyr.append(position)
                elif params['type'] == 'topological domain':
                    # The current feture contains topological domain information (information about the subcellular compartment where each non-membrane region of a membrane-spanning protein is found).
                    domain = params['description']  # Record the domain.
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein this domain occurs.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the domain, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the domain, then record this.
                        endLoc = ''
                    topologicalDomain.append(domain + ',' + beginLoc + ',' + endLoc)
                elif params['type'] == 'signal peptide':
                    # The current feature contains information about a signal peptide.
                    if params.has_key('evidence'):
                        # The signal peptide location has been experimentally proven.
                        evidence = 'Y'
                    else:
                        # The signal peptide location is not recorded as having been experimentally proven.
                        evidence = 'N'
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein the signal peptide is found.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the signal peptide, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the signal peptide, then record this.
                        endLoc = ''
                    signalPeptide.append(beginLoc + ',' + endLoc + ',' + evidence)
                elif (params['type'] == 'transmembrane region' and params.has_key('description') and params['description'][:7] == 'Helical'):
                    # The current feature contains information about an alpha helical transmembrane region.
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein the transmembrane region is found.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the transmembrane region, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the transmembrane region, then record this.
                        endLoc = ''
                    transMem.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'turn':
                    # The current feature contains information about a hydrogen bonded turn.
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein the turn is found.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the turn, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the turn, then record this.
                        endLoc = ''
                    turns.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'helix':
                    # The current feature contains information about an alpha helix.
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein the alpha helix is found.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the alpha helix, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the alpha helix, then record this.
                        endLoc = ''
                    helices.append(beginLoc + ',' + endLoc)
                elif params['type'] == 'strand':
                    # The current feature contains information about a beta strand.
                    locationTag = j.children[0]  # Get information about between which two amino acids in the protein the beta strand is found.
                    try:
                        beginLoc = locationTag.children[0].parameters['position']
                    except:
                        # If there's no information about the start site of the beta strand, then record this.
                        beginLoc = ''
                    try:
                        endLoc = locationTag.children[1].parameters['position']
                    except:
                        # If there's no information about the end site of the beta strand, then record this.
                        endLoc = ''
                    strands.append(beginLoc + ',' + endLoc)

            # Add the feature information to the tuple of information about the protein.
            proteinRecord[38] = ';'.join(set(glycosylationO)) if glycosylationO != [] else 'NA'
            proteinRecord[39] = ';'.join(set(glycosylationN)) if glycosylationN != [] else 'NA'
            proteinRecord[40] = ';'.join(set(phosphSer)) if phosphSer != [] else 'NA'
            proteinRecord[41] = ';'.join(set(phosphThr)) if phosphThr != [] else 'NA'
            proteinRecord[42] = ';'.join(set(phosphTyr)) if phosphTyr != [] else 'NA'
            proteinRecord[44] = ';'.join(set(topologicalDomain)) if topologicalDomain != [] else 'NA'
            proteinRecord[46] = ';'.join(set(signalPeptide)) if signalPeptide != [] else 'NA'
            proteinRecord[47] = ';'.join(set(transMem)) if transMem != [] else 'NA'
            proteinRecord[48] = ';'.join(set(turns)) if turns != [] else 'NA'
            proteinRecord[49] = ';'.join(set(helices)) if helices != [] else 'NA'
            proteinRecord[50] = ';'.join(set(strands)) if strands != [] else 'NA'

            # Add the information about the protein sequence to the record tuple.
            sequence = (tags[i]['uniprot.entry.sequence'][0].data).replace(' ', '')
            proteinRecord[-1] = sequence

            # Add the protein's record tuple to the list of protein record tuples.
            proteinTuples.append(tuple(proteinRecord))

    # Write out the files generated from the UniProt data.
    utilities.list2file.main(UPRepresentativeAccessions, UPHumanAccessions)
    utilities.list2file.main(UPAccessionMapping, UPHumanAccessionMap)
    utilities.list2file.main(UPProteinNames, UPHumanNames)
    utilities.list2file.main(UPDrugMapping, UPDrugIDs)
    utilities.list2file.main(proteinTuples, UPProteinInfo)
    utilities.list2file.main(proteinExternalLinks, UPExternalLinks)
    utilities.list2file.main(PPIData, UPPPIData)

def split_uniprot(XMLInputFile, entriesPerSplit, folderUP):
    """
    Splits the UniProt XML file into smaller XML files containing entriesPerSplit UniProt entries per file.
    """

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