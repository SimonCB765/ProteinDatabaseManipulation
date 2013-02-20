'''
Created on 14 Oct 2011

@author: Simon Bull
'''

import utilities.list2file

def main(ensemblTranscripts, ensemblParsedTranscripts, ensemblGermSNPResults, ensemblParsedGermVariants):

    parse_variants(ensemblGermSNPResults, ensemblParsedGermVariants)
    parse_transcript(ensemblTranscripts, ensemblParsedTranscripts)

def parse_variants(ensemblGermSNPResults, ensemblParsedGermVariants):
    """
    Takes a file of germline mutations and the genes that they occur in, and returns a file of mutations and their types.
    ensemblParsedGermVariants - A file of 25-tuples, with one on each line.
        The first element is the Ensembl transcript ID.
        The second element is mutation ID.
        The third element is the Ensembl gene ID that the transcript in the first element comes from.
        The fourth element is the amino acid change cause by the mutation ('NA' if no change is caused/known).
        Elements five through twenty-five indicate the consequence of the mutation identified in the second element on the transcript in the first element. A 1 for an element indicates that the consequence corresponding to the element occurs.
            The order of the primary sites in the tuple is the same as the dictionary consequenceDict in the code.
    """

    # The possible types of mutation, and their index in the mutation information tuple.
    consequenceDict = {'3PRIME_UTR' : 4, '5PRIME_UTR' : 5, 'CODING_UNKNOWN' : 6, 'COMPLEX_INDEL' : 7,
                       'DOWNSTREAM' : 8, 'ESSENTIAL_SPLICE_SITE' : 9, 'FRAMESHIFT_CODING' : 10, 'INTERGENIC' : 11,
                       'INTRONIC' : 12, 'NMD_TRANSCRIPT' : 13, 'NON_SYNONYMOUS_CODING' : 14, 'PARTIAL_CODON' : 15,
                       'REGULATORY_REGION' : 16, 'SPLICE_SITE' : 17, 'STOP_GAINED' : 18, 'STOP_LOST' : 19,
                       'SYNONYMOUS_CODING' : 20, 'TRANSCRIPTION_FACTOR_BINDING_MOTIF' : 21, 'UPSTREAM' : 22,
                       'WITHIN_MATURE_miRNA' : 23, 'WITHIN_NON_CODING_GENE' : 24}
    # The initialised empty tuple.
    defaultTuple = ['trans', 'variant', 'gene', 'NA', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    readIn = open(ensemblGermSNPResults, 'r')
    writeOut = open(ensemblParsedGermVariants, 'w')
    currentGene = ''
    geneSNPTransPairs = {}  # Records the tuples of the mutation information for each gene.
    for line in readIn:
        chunks = line[:-1].split('\t')
        geneID = chunks[0]  # The Ensembl gene ID of the gene involved in this mtuation relation is the first element of the line.
        if geneID != currentGene:
            # A new gene has been found (the file is ordered by genes, so all mutation related to one gene occur on contiguous lines).
            tuplesToWrite = set([geneSNPTransPairs[i]['tuple'] for i in geneSNPTransPairs.keys() if geneSNPTransPairs[i]['valid']])  # Generate the tuple to write out for the last gene encountered.
            for i in tuplesToWrite:
                # Write out all the tuples for the last gene encountered.
                writeOut.write(str(i) + '\n')
            geneSNPTransPairs = {}
            currentGene = geneID

        # The gene on the current line is the same as the one on the last line.
        currentTuple = list(defaultTuple)
        currentTuple[0] = chunks[1]  # Set the Ensembl transcript ID in the current tuple.
        currentTuple[1] = chunks[2]  # Set the ID of the mutation found on the line.
        currentTuple[2] = geneID  # Set the Ensembl gene ID of the gene found on the line.
        changeInAA = chunks[3]  # Record the change in the amino acids in the protein sequence caused by the mutation.
        currentTuple[3] = 'NA' if changeInAA == '' else changeInAA  # If no change in amino acids was recorded on the line, then set the change to 'NA'.
        consequence = chunks[4].split(',')  # There can be multiple consequences of the mutation (UPSTREAM, DOWNSTREAM, ESSENTIAL_SPLICE_SITE, etc.), all separated by commas.
        for j in consequence:
            # Set all entries in the tuple for the consequences observed for the mutation to 1. All other consequences remain 0.
            currentTuple[consequenceDict[j]] = 1

        # The unique key for each mutation is the mutation ID combined with the transcript ID. This is because a mutation can affect multiple transcripts, and therefore will not be unique.
        key = tuple([transcriptID, variantID])
        if not geneSNPTransPairs.has_key(key):
            # If the mutation ID transcript ID pair has not been seen before, then reocrd it for the first time.
            geneSNPTransPairs[key] = {'valid' : True, 'cons' : consequence, 'tuple' : tuple(currentTuple)}
        else:
            if geneSNPTransPairs[key]['valid'] == False:
                # If the (SNP, transcript) pair has already been marked as invalid, then ignore it this time.
                continue
            elif geneSNPTransPairs[key]['cons'] != consequence:
                # There are two (or more) occurunces of the (SNP, transcript) pair, but the consequence(s) are not the same.
                # Therefore, the mutation ID transcript ID pair should be marked as invalid.
                geneSNPTransPairs[key]['valid'] = False
    writeOut.close()
    readIn.close()

def parse_transcript(ensemblTranscripts, ensemblParsedTranscripts):
    """
    Takes a file containing genes, the number of transcripts they have and the type of those transcripts. Returns a file of gene transcript information.
    ensemblParsedTranscripts - A file of 6-tuples, with one on each line.
        The first element is the Ensembl gene ID.
        The second element is the number of different transcripts that the gene in the first element generates.
        The third element is the number of transcripts the gene in the first element generates that are protein coding.
        The fourth element is the number of transcripts the gene in the first element generates that are retain intron.
        The fifth element is the number of transcripts the gene in the first element generates that are processed transcript.
        The sixth element is the number of transcripts the gene in the first element generates that are nonsense mediated decay.
    """

    parsedOutput = {}
    readIn = open(ensemblTranscripts, 'r')
    for line in readIn:
        chunks = line[:-1].split('\t')
        ensemblGeneID = chunks[0]
        transcriptCount = int(chunks[2])
        codesForProtein = 1 if chunks[3] == 'protein_coding' else 0
        retainedIntron = 1 if chunks[3] == 'retained_intron' else 0
        processedTranscript = 1 if chunks[3] == 'processed_transcript' else 0
        nonsenseMediatedDecay = 1 if chunks[3] == 'nonsense_mediated_decay' else 0
        if parsedOutput.has_key(ensemblGeneID):
            parsedOutput[ensemblGeneID]['ProteinCoding'] += codesForProtein
            parsedOutput[ensemblGeneID]['RetainedIntron'] += retainedIntron
            parsedOutput[ensemblGeneID]['ProcessedTranscript'] += processedTranscript
            parsedOutput[ensemblGeneID]['NonsenseMediatedDecay'] += nonsenseMediatedDecay
        else:
            parsedOutput[ensemblGeneID] = {}
            parsedOutput[ensemblGeneID]['Count'] = transcriptCount
            parsedOutput[ensemblGeneID]['ProteinCoding'] = codesForProtein
            parsedOutput[ensemblGeneID]['RetainedIntron'] = retainedIntron
            parsedOutput[ensemblGeneID]['ProcessedTranscript'] = processedTranscript
            parsedOutput[ensemblGeneID]['NonsenseMediatedDecay'] = nonsenseMediatedDecay
    readIn.close()

    parsedOutput = [tuple([i, parsedOutput[i]['Count'], parsedOutput[i]['ProteinCoding'], parsedOutput[i]['RetainedIntron'],
                           parsedOutput[i]['ProcessedTranscript'], parsedOutput[i]['NonsenseMediatedDecay']])
                    for i in parsedOutput.keys()]

    utilities.list2file.main(parsedOutput, ensemblParsedTranscripts)

def update_xref_and_ensembl_IDs(ensemblExternalIDsOne, ensemblExternalIDsTwo, uniprotExternalIDs, ensemblGeneFile):
    """
    Takes two files containing the cross-referencing of representative human UniProt accessions with Ensembl gene IDs, Ensembl transcript IDs, Ensembl protein IDs, Entrez Gene IDs, UniGene cluster IDs and HGNC gene IDs.
    Returns two files.
        uniprotExternalIDs contains the cross-referencing of the represetnative UniProt human accessions with external databases.
        ensemblGeneFile contains all the Ensembl gene IDs for Ensembl genes that are linked to a representative UniProt human accession.
    uniprotExternalIDs - A comma separated (csv) file, with six elements on each line.
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
    ensemblGeneFile - A file with one Ensembl gene ID on each line.
    """

    # Record the new cross-reference information extracted from Ensembl.
    ensemblXref = {}
    readEnsembl = open(ensemblExternalIDsOne, 'r')
    for line in readEnsembl:
        line = line.strip()
        chunks = line.split('\t')
        EnsemblID = chunks[0]  # The first entry on a line is the Ensembl gene ID.
        EnsemblTranscript = chunks[1]  # The second entry on a line is the Ensembl transcript ID.
        EnsemblProtein = chunks[2]  # The third entry on a line is the Ensembl protein ID.
        EnsemblEntrezGeneID = chunks[3]  # The fourth entry on a line is the Entrez Gene ID.
        EnsemblUnigeneID = chunks[4]  # The fifth entry on a line is the UniGene cluster ID.
        EnsemblUPAcc = chunks[5]  # The sixth entry on a line is the UniProt accession.

        if not ensemblXref.has_key(EnsemblUPAcc):
            # If the UniProt accession has not already been recorded, then add an empty record for it.
            # Each record is a set to prevent duplicate cross-references being recorded.
            ensemblXref[EnsemblUPAcc] = {'EnsemblGenes' : set([]), 'Gene' : set([]), 'UniGene' : set([]), 'HGNC' : set([]),
                                         'EnsemblTriplets' : set([])}

        ensemblXref[EnsemblUPAcc]['EnsemblGenes'].add(EnsemblID)  # Add the Ensembl gene ID to the set of Ensembl gene IDs that are linked to the UniProt accession.
        ensemblXref[EnsemblUPAcc]['EnsemblTriplets'].add(EnsemblID + '-' + EnsemblTranscript + '-' + EnsemblProtein)  # Add the Ensembl triplet of IDs to the set of Ensembl triplets of IDs that are linked to the UniProt accession.
        ensemblXref[EnsemblUPAcc]['Gene'].add(EnsemblEntrezGeneID)  # Add the Entrez Gene ID to the set of Entrez Gene IDs that are linked to the UniProt accession.
        ensemblXref[EnsemblUPAcc]['UniGene'].add(EnsemblUnigeneID)  # Add the UniGene cluster ID to the set of UniGene cluster IDs that are linked to the UniProt accession.
    readEnsembl.close()
    readEnsembl = open(ensemblExternalIDsTwo, 'r')
    for line in readEnsembl:
        line = line.strip()
        chunks = line.split('\t')
        if len(chunks) < 5:
            # As the only new information here is the HGNC gene ID, we don;t care about the information on the line if the HGNC gene ID is not present.
            # If there are less than 5 splits the HGNC gene ID information is absent, and therefore the line can be skipped.
            continue
        EnsemblID = chunks[0]  # The first entry on a line is the Ensembl gene ID.
        EnsemblTranscript = chunks[1]  # The second entry on a line is the Ensembl transcript ID.
        EnsemblProtein = chunks[2]  # The third entry on a line is the Ensembl protein ID.
        EnsemblUPAcc = chunks[3]  # The fourth entry on a line is the UniProt accession.
        EnsemblHGNCID = chunks[4]  # The fifth entry on a line is the HGNC gene ID.

        if not ensemblXref.has_key(EnsemblUPAcc):
            # If the UniProt accession has not already been recorded, then add an empty record for it.
            # Each record is a set to prevent duplicate cross-references being recorded.
            ensemblXref[EnsemblUPAcc] = {'EnsemblGenes' : set([]), 'Gene' : set([]), 'UniGene' : set([]), 'HGNC' : set([]),
                                         'EnsemblTriplets' : set([])}

        ensemblXref[EnsemblUPAcc]['EnsemblGenes'].add(EnsemblID)  # Add the Ensembl gene ID to the set of Ensembl gene IDs that are linked to the UniProt accession.
        ensemblXref[EnsemblUPAcc]['EnsemblTriplets'].add(EnsemblID + '-' + EnsemblTranscript + '-' + EnsemblProtein)  # Add the Ensembl triplet of IDs to the set of Ensembl triplets of IDs that are linked to the UniProt accession.
        ensemblXref[EnsemblUPAcc]['HGNC'].add(EnsemblHGNCID)  # Add the HGNC gene ID to the set of HGNC gene IDs that are linked to the UniProt accession.
    readEnsembl.close()

    # Record the cross-reference information extracted from UniProt.
    uniprotXref = {}
    readUP = open(uniprotExternalIDs, 'r')
    for line in readUP:
        line = line.strip()
        chunks = line.split(',')
        UPAcc = chunks[0]
        UPGeneID = chunks[1].split(';')
        UPUnigeneID = chunks[2].split(';')
        UPGOID = chunks[3].split(';')
        UPHGNCID = chunks[4].split(';')
        uniprotXref[UPAcc] = {'Gene' : set(UPGeneID), 'UniGene' : set(UPUnigeneID), 'GO' : set(UPGOID), 'HGNC' : set(UPHGNCID)}
    readUP.close()

    ensemblGeneIDs = set([])  # Record the set of Ensembl gene IDs that are linked to the representative human UniProt accessions.
    newXrefs = []
    for i in ensemblXref.keys():
        # Remove the empty string if it is present in any of the entries.
        ensemblXref[i]['Gene'] -= set([''])
        uniprotXref[i]['Gene'] -= set([''])
        ensemblXref[i]['UniGene'] -= set([''])
        uniprotXref[i]['UniGene'] -= set([''])
        uniprotXref[i]['GO'] -= set([''])
        ensemblXref[i]['HGNC'] -= set([''])
        uniprotXref[i]['HGNC'] -= set([''])
        ensemblXref[i]['EnsemblGenes'] -= set([''])
        ensemblXref[i]['EnsemblTriplets'] -= set([''])

        # Add all the Ensembl gene IDs to the set of Ensembl gene IDs that are linked to the representative human UniProt accessions.
        ensemblGeneIDs |= ensemblXref[i]['EnsemblGenes']

        # For each protein record the union of all the external database cross-references extracted from both Ensembl and UniProt.
        output = [i]
        output.append(';'.join(ensemblXref[i]['Gene'].union(uniprotXref[i]['Gene'])))
        output.append(';'.join(ensemblXref[i]['UniGene'].union(uniprotXref[i]['UniGene'])))
        output.append(';'.join(uniprotXref[i]['GO']))
        output.append(';'.join(ensemblXref[i]['HGNC'].union(uniprotXref[i]['HGNC'])))
        output.append(';'.join(ensemblXref[i]['EnsemblTriplets']))
        output = ','.join(output)
        newXrefs.append(output)

    utilities.list2file.main(newXrefs, uniprotExternalIDs)
    utilities.list2file.main(list(ensemblGeneIDs), ensemblGeneFile)
