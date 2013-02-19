import re

def main(COSMICData, COSMICParsedGene, COSMICParsedMutation, COSMICParsedGene2Mutation):
    """
    Takes a file containing the COSMIC data, and returns three files.
        COSMICParsedGene contains a mapping of genes to the sites where they have been observed as being expressed.
        COSMICParsedMutation ccontains mutations along with the amino acid changes they cause, the type of the mutation and the sites where they have been found.
        COSMICParsedGene2Mutation contains a mapping of genes to mutations that occur in them.
    COSMICParsedGene - A file of tuples with forty-nine elements, with on tuple on each line.
        The first element is the gene ID.
        The second element is the Ensembl transcript ID that corresponds to the gene ID in the first element.
        The third element is the HGNC gene ID that corresponds to the gene ID in the first element.
        Elements four through forty-nine indicate the number of times that the gene in the first element was observed in each of the forty-six primary sites.
            The order of the primary sites in the tuple is the same as the list primarySites in the code.
    COSMICParsedMutation - A file of tuples with sixty-five elements, with on tuple on each line.
        The first element is the mutation ID.
        The second element is the amino acid change.
        The third element is 'Somatic', 'Unknown' or 'Germline' depending on whether the mutation is somatic, unknown or germline respectively.
        Elements four through twenty indicate the type of the mutation in the first element. Only one of the sixteen options will be set to 1 (indicating the mutation is of that type), the rest will be set to 0.
            The order of the mutation types in the tuple is the same as the list mutationTypes in the code.
        Elements twenty-one through sixty-five indicate the number of times that the mutation in the first element was observed in each of the forty-six primary sites.
            The order of the primary sites in the tuple is the same as the list primarySites in the code.
    COSMICParsedGene2Mutation - A file of tuples with two elements, with one tuple on each line.
        The first element is the gene ID.
        The second element is the mutation ID.
    """

    # The sites where mutation can occur.
    primarySites = ['adrenal_gland', 'autonomic_ganglia', 'biliary_tract', 'bone', 'breast', 'central_nervous_system',
                    'cervix', 'endometrium', 'eye', 'fallopian_tube', 'female_genital_tract_(site_indeterminate)',
                    'gastrointestinal_tract_(site_indeterminate)', 'genital_tract', 'haematopoietic_and_lymphoid_tissue',
                    'kidney', 'large_intestine', 'liver', 'lung', 'mediastinum', 'meninges', 'midline_organs',
                    'oesophagus', 'ovary', 'pancreas', 'paratesticular_tissues', 'parathyroid', 'penis', 'pericardium',
                    'peritoneum', 'pituitary', 'placenta', 'pleura', 'prostate', 'retroperitoneum', 'salivary_gland',
                    'skin', 'small_intestine', 'soft_tissue', 'stomach', 'testis', 'thymus', 'thyroid',
                    'upper_aerodigestive_tract', 'urinary_tract', 'vagina', 'vulva'
                    ]

    # The types of mutation that can occur.
    mutationTypes = ['Complex', 'Complex - compound substitution', 'Complex - deletion inframe', 'Complex - frameshift', 'Complex - insertion inframe',
                     'Deletion - In frame', 'Deletion - Frameshift',
                     'Insertion - Frameshift', 'Insertion - In frame',
                     'No detectable mRNA/protein', 'Nonstop extension',
                     'Substitution - coding silent', 'Substitution - Missense', 'Substitution - Nonsense',
                     'Unknown',
                     'Whole gene deletion'
                     ]

    # The three different possibilities for classifying mutations.
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
            # If there is an Ensembl transcript ID in the accession column.
            acc = match.group(0)
        else:
            # No accession of interest to record.
            acc = ''
        hgnc = chunks[2]
        if hgnc != '':
            # If there is an HGNC gene ID recorded for the gene on the line, then record the HGNC gene ID.
            hgnc = int(hgnc)
        else:
            hgnc = -1
        primarySite = chunks[6]  # Record the primary mutation site.
        mutationID = chunks[11]  # Record the ID of the mutation.
        aaChange = chunks[13]  # Record the change in the amino acid sequence caused by the mutation.
        mutationType = chunks[14]  # Record the type of the mutation.
        somatic = chunks[20]  # Record whether the mutation is somatic, germline or unknown.
        if somatic in somaticMutation:
            somatic = 'Somatic'
        elif somatic in germlineMutation:
            somatic = 'Germline'
        else:
            somatic = 'Unknown'

        if geneData.has_key(gene):
            # If the gene has already been encountered.
            if primarySite in primarySites:
                # If the primary site is in one of the sites of interest, then increment the number of times that this gene was seen at the site.
                geneData[gene]['Site'][primarySites.index(primarySite)] += 1
        else:
            # If the gene has not been encountered before, then initialise all the sites to 0 occurences.
            # Set the primary site encountered on the current line to 1.
            startPrimarySite = [0 if i != primarySite else 1 for i in primarySites]
            geneData[gene] = {'Acc' : acc, 'HGNC' : hgnc, 'Site' : startPrimarySite}

        if mutationID != '':
            # If there is a mutation.
            mutationID = int(mutationID)
            gene2Mut.add(tuple([gene, mutationID]))
            if mutationData.has_key(mutationID):
                # If the mutation has been seen before.
                if primarySite in primarySites:
                    # If the primary site is in one of the sites of interest, then increment the number of times that this mutation was seen at the site.
                    mutationData[mutationID]['Site'][primarySites.index(primarySite)] += 1
            else:
                # If the mutation has not been encountered before, then initialise all the sites to 0 occurences.
                # Set the primary site encountered on the current line to 1.
                # also indicate the type of mutation by setting it to 1 and all other mutation types to 0.
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