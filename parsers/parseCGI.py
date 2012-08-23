import re

def main(CGIData, CGIHGNCParsed, CGIUPAccParsed):

    readIn = open(CGIData, 'r')
    validGene = False
    newSentence = False
    humanGene = False
    validSentence = False
    evidenceExists = False
    tempHGNC = ''
    tempUPAcc = ''
    cancerHGNCs = set([])
    cancerUPAccs = set([])
    for line in readIn:
        if '<GeneEntry>' in line:
            if validGene:
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
                        humanGene = True
                elif matchValidSentence:
                    validSentence = True
                elif matchEvidence:
                    evidence = matchEvidence.group(0)
                    if not evidence in ['false positive', 'not_assigned', 'EV-AS-NAS', 'EV-COMP', 'EV-COMP-AINF']:
                        evidenceExists = True
        elif validGene:
            continue
        else:
            matchHGNC = re.search('(?<=<HgncID>).*(?=</HgncID>)', line)
            matchUP = re.search('(?<=<UniProtID>).*(?=</UniProtID>)', line)
            if matchHGNC:
                tempHGNC = int(matchHGNC.group(0))
            elif matchUP:
                tempUPAcc = matchUP.group(0)
            elif re.search('<Sentence>', line):
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