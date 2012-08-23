import re

def main(COSMICData, COSMICParsedGene, COSMICParsedMutation, COSMICParsedGene2Mutation):

    primarySites = ['adrenal_gland', 'autonomic_ganglia', 'biliary_tract', 'bone', 'breast', 'central_nervous_system',
                    'cervix', 'endometrium', 'eye', 'fallopian_tube', 'female_genital_tract_(site_indeterminate)',
                    'gastrointestinal_tract_(site_indeterminate)', 'genital_tract', 'haematopoietic_and_lymphoid_tissue',
                    'kidney', 'large_intestine', 'liver', 'lung', 'mediastinum', 'meninges', 'midline_organs',
                    'oesophagus', 'ovary', 'pancreas', 'paratesticular_tissues', 'parathyroid', 'penis', 'pericardium',
                    'peritoneum', 'pituitary', 'placenta', 'pleura', 'prostate', 'retroperitoneum', 'salivary_gland',
                    'skin', 'small_intestine', 'soft_tissue', 'stomach', 'testis', 'thymus', 'thyroid',
                    'upper_aerodigestive_tract', 'urinary_tract', 'vagina', 'vulva'
                    ]
    mutationTypes = ['Complex', 'Complex - compound substitution', 'Complex - deletion inframe', 'Complex - frameshift', 'Complex - insertion inframe',
                     'Deletion - In frame', 'Deletion - Frameshift',
                     'Insertion - Frameshift', 'Insertion - In frame',
                     'No detectable mRNA/protein', 'Nonstop extension',
                     'Substitution - coding silent', 'Substitution - Missense', 'Substitution - Nonsense',
                     'Unknown',
                     'Whole gene deletion'
                     ]
    somaticMutation = ['Confirmed somatic variant', 'Reported in another cancer sample as somatic']
    germlineMutation = ['Confirmed germline variant', 'Reported in another sample as germline']
    unknownMutation = ['Not specified', 'Variant of unknown origin']

    geneData = {}
    gene2Mut = set([])
    mutationData = {}
    readIn = open(COSMICData, 'r')
    header = readIn.readline()  # Remove the header line.
    for line in readIn:
        chunks = line.split('\t')
        gene = chunks[0]
        acc = chunks[1]
        match = re.search('ENST[0-9]{11}', acc)
        if match:
            acc = match.group(0)
        else:
            acc = ''            
        hgnc = chunks[2]
        if hgnc != '':
            hgnc = int(hgnc)
        else:
            hgnc = -1
        primarySite = chunks[6]
        mutationID = chunks[11]
        aaChange = chunks[13]
        mutationType = chunks[14]
        somatic = chunks[20]
        if somatic in somaticMutation:
            somatic = 'Somatic'
        elif somatic in germlineMutation:
            somatic = 'Germline'
        else:
            somatic = 'Unknown'

        if geneData.has_key(gene):
            if primarySite in primarySites:
                geneData[gene]['Site'][primarySites.index(primarySite)] += 1
        else:
            startPrimarySite = [0 if i != primarySite else 1 for i in primarySites]
            geneData[gene] = {'Acc' : acc, 'HGNC' : hgnc, 'Site' : startPrimarySite}

        if mutationID != '':
            # If there is a mutation.
            mutationID = int(mutationID)
            gene2Mut.add(tuple([gene, mutationID]))

            if mutationData.has_key(mutationID):
                if primarySite in primarySites:
                    mutationData[mutationID]['Site'][primarySites.index(primarySite)] += 1
            else:
                startPrimarySite = [0 if i != primarySite else 1 for i in primarySites]
                startMutationType = [0 if i != mutationType else 1 for i in mutationTypes]
                mutationData[mutationID] = {'AA' : aaChange, 'Somatic' : somatic, 'MutType' : startMutationType, 'Site' : startPrimarySite}
    readIn.close()

    writeOut = open(COSMICParsedGene, 'w')
    for i in geneData.keys():
        writeOut.write(str(tuple([i, geneData[i]['Acc'], geneData[i]['HGNC']] + geneData[i]['Site'])) + '\n')
    writeOut.close()

    writeOut = open(COSMICParsedMutation, 'w')
    for i in mutationData.keys():
        writeOut.write(str(tuple([i, mutationData[i]['AA'], mutationData[i]['Somatic']] + mutationData[i]['MutType'] + mutationData[i]['Site'])) + '\n')
    writeOut.close()

    writeOut = open(COSMICParsedGene2Mutation, 'w')
    for i in gene2Mut:
        writeOut.write(str(i) + '\n')
    writeOut.close()