'''
Created on 8 Nov 2011

@author: Simonial
'''

import subprocess

import utilities.list2file

def main(ensemblHomologScript, ensemblTaxonomyMap, ensemblHomologData, ensemblGenomesHomologScript,
         ensemblGenomesTaxonomyMap, ensemblGenomesHomologData, ensemblParsedHomology, ensemblGeneIDs):
    """
    Takes files containing the Ensembl Perl API scripts and the list of Ensembl gene IDs linked to UniProt represetnative accessions (along with temporary files for the storing of the homolg data).
    Returns a files containing the homolog information for all the Ensembl gene IDs linked to a representative UniProt accession, and two files containing the mappings of the taxonomy IDs to organisms used in the Esnembl databases.
    ensemblTaxonomyMap - A tab separated (tsv) file, with two elements on each line.
        The first element is the organism ID.
        The second element is the name of the organism that corresponds to the ID in the first element.
    ensemblGenomesTaxonomyMap - A tab separated (tsv) file, with two elements on each line.
        The first element is the organism ID.
        The second element is the name of the organism that corresponds to the ID in the first element.
    ensemblParsedHomology - A file containing 10-tuples, with one on each line.
        The first element is the human Ensembl gene ID.
        The second element is the ID of the gene that is homologus to the gene in the first element.
        The third element is the organism that the gene in the second element comes from.
        The fourth element is the type of homolog.
        The fifth element is the most recent common ancestor of the genes in the first two elements.
        The sixth element is a value for dN (non-synonymous substitutions per non-synonymous site) (-1.0 if a value is not present in Ensembl for this).
        The seventh element is a value for dS (synonymous substitutions per synonymous site) (-1.0 if a value is not present in Ensembl for this).
        The eighth element is the percentage of the peptide which has been aligned (-1.0 if a value is not present in Ensembl for this).
        The ninth element is the percentage of identity between both homologs (-1.0 if a value is not present in Ensembl for this).
        The tenth element is the percentage of positivity (similarity) between both homologs (-1.0 if a value is not present in Ensembl for this).
    """

    # Extract the Ensembl homolog information.
    subprocess.call(['perl', ensemblHomologScript, ensemblGeneIDs, ensemblHomologData, ensemblTaxonomyMap])

    # Extract the Ensembl Genomes pan-taxonomic homologs.
    subprocess.call(['perl', ensemblGenomesHomologScript, ensemblGeneIDs, ensemblGenomesHomologData, ensemblGenomesTaxonomyMap])

    # Parse the information from the Ensembl database (the multi species data).
    ensemblTaxonomyDict = {}
    readTaxon = open(ensemblTaxonomyMap, 'r')
    for line in readTaxon:
        chunks = (line.rstrip()).split('\t')
        ensemblTaxonomyDict[chunks[0]] = chunks[1]
    readTaxon.close()

    ensemblHomologDict = {}
    readHomologs = open(ensemblHomologData, 'r')
    for line in readHomologs:
        chunks = (line.rstrip()).split('\t')
        ensemblHomologDict[tuple([chunks[0], chunks[2]])] = tuple([chunks[0], chunks[2], ensemblTaxonomyDict[chunks[3]],
                                                                   chunks[4], chunks[5], float(chunks[6]),
                                                                   float(chunks[7]), float(chunks[8]),
                                                                   float(chunks[9]), float(chunks[10])])
    readHomologs.close()

    # Parse the information from the Ensembl Genome database (the pan_homology data).
    ensemblGenomeTaxonomyDict = {}
    readTaxon = open(ensemblGenomesTaxonomyMap, 'r')
    for line in readTaxon:
        chunks = (line.rstrip()).split('\t')
        ensemblGenomeTaxonomyDict[chunks[0]] = chunks[1]
    readTaxon.close()

    ensemblGenomeHomologDict = {}
    readHomologs = open(ensemblGenomesHomologData, 'r')
    for line in readHomologs:
        chunks = (line.rstrip()).split('\t')
        ensemblGenomeHomologDict[tuple([chunks[0], chunks[2]])] = tuple([chunks[0], chunks[2],
                                                                         ensemblGenomeTaxonomyDict[chunks[3]],
                                                                         chunks[4], chunks[5], float(chunks[6]),
                                                                         float(chunks[7]), float(chunks[8]),
                                                                         float(chunks[9]), float(chunks[10])])
    readHomologs.close()

    # Output the parse homology data. Preference is given to data extracted from Ensembl Genome as it is (potentially) more up to date.
    homologs = []
    ensemblGenomeHomologs = set(ensemblGenomeHomologDict.keys())
    for i in ensemblGenomeHomologs:
        homologs.append(ensemblGenomeHomologDict[i])
    for i in ensemblHomologDict.keys():
        if i in ensemblGenomeHomologs:
            # If the homolog was found in the Ensembl Genomes data, then ifgnore it here.
            continue
        else:
            homologs.append(ensemblHomologDict[i])

    utilities.list2file.main(homologs, ensemblParsedHomology)