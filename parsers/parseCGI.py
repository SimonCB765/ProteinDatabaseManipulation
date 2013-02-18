import re

def main(CGIData, CGIHGNCParsed, CGIUPAccParsed):
    """
    Takes the file containg the cancer gene index data, and returns one file containing UniProt accessions and one containing HGNC gene IDs.
    CGIHGNCIDs (CGIHGNCIDs) - A file with one HGNC gene ID on each line.
    CGIUPAccessions (CGIUPAccessions) - A file with one UniProt accession on each line.
    """

    readIn = open(CGIData, 'r')
    validGene = False  # Records whether the gene found is valid for our purposes. for this to be true the gene must be human, there must be evidence implicating it in cancer and the <sentence> tag that provides the information must be 'finished'.
    newSentence = False  # Whether we have foudn the start of a <sentence> tag.
    humanGene = False  # Whether the current gene is human.
    validSentence = False  # Whether the sentence being examined is 'finished'.
    evidenceExists = False  # Whether evidence exists for the current gene's implication in cancer.
    tempHGNC = ''
    tempUPAcc = ''
    cancerHGNCs = set([])  # Records all the HGNC IDs that are implicated in cancer.
    cancerUPAccs = set([])  # Records all the UniProt accessions that are implicated in cancer.
    for line in readIn:
        if '<GeneEntry>' in line:
            # If a new gene entry has been found.
            if validGene:
                # If the last gene was a valid one, then add the HGNC gene ID and UniProt accession to the list of those implicated in cancer.
                cancerHGNCs.add(tempHGNC)
                cancerUPAccs.add(tempUPAcc)
            validGene = False
            newSentence = False
            humanGene = False
            validSentence = False
            evidenceExists = False
            tempHGNC = ''
            tempUPAcc = ''
        elif newSentence:
            if re.search('</Sentence>', line):
                newSentence = False
                if validSentence and humanGene and evidenceExists:
                    validGene = True
            else:
                matchOrganism = re.search('(?<=<Organism>).*(?=</Organism>)', line)
                matchValidSentence = re.search('<SentenceStatusFlag>finished</SentenceStatusFlag>', line)
                matchEvidence = re.search('(?<=<EvidenceCode>).*(?=</EvidenceCode>)', line)
                if matchOrganism:
                    if matchOrganism.group(0).lower() == 'human':
                        # The organism that the gene occurs in is human.
                        humanGene = True
                elif matchValidSentence:
                    # The <sentence> tag being investigated is 'finished', and therefore contains useable information.
                    validSentence = True
                elif matchEvidence:
                    evidence = matchEvidence.group(0)
                    if not evidence in ['', 'false positive', 'not_assigned', 'EV-AS-NAS', 'EV-COMP'] or evidence[:12] == 'EV-COMP-AINF':
                        # If the evidence code is not 'false positive', 'not_assigned', 'EV-AS-NAS' or 'EV-COMP', is not empty and does not begin with 'EV-COMP-AINF', then
                        # the evidence for the implication of the gene/protein in cancer is deemed to be acceptable.
                        evidenceExists = True
        elif validGene:
            # If you have found that the current gene is validly implicated in cancer, then cycle through to the end of the gene entry.
            continue
        else:
            matchHGNC = re.search('(?<=<HgncID>).*(?=</HgncID>)', line)
            matchUP = re.search('(?<=<UniProtID>).*(?=</UniProtID>)', line)
            if matchHGNC:
                # The line contains an HGNC gene ID.
                tempHGNC = int(matchHGNC.group(0))
            elif matchUP:
                # The line contains a UniProt accession.
                tempUPAcc = matchUP.group(0)
            elif re.search('<Sentence>', line):
                # The start of a sentence tag has been found.
                newSentence = True
    readIn.close()
    cancerHGNCs -= set([''])
    cancerUPAccs -= set([''])

    writeTo = open(CGIHGNCParsed, 'w')
    for i in cancerHGNCs:
        writeTo.write(str(i) + '\n')
    writeTo.close()

    writeTo = open(CGIUPAccParsed, 'w')
    for i in cancerUPAccs:
        writeTo.write(i + '\n')
    writeTo.close()