'''
Created on 8 Nov 2011

@author: Simonial
'''

import subprocess

import utilities.list2file

def main(ensemblHomologScript, ensemblTaxonomyMap, ensemblHomologData, ensemblGenomesHomologScript,
         ensemblGenomesTaxonomyMap, ensemblGenomesHomologData, ensemblParsedHomology, ensemblGeneIDs):

##    # Extract the Ensembl homolog information.
##    subprocess.call(['perl', ensemblHomologScript, ensemblGeneIDs, ensemblHomologData, ensemblTaxonomyMap])
    
##    # Extract the Ensembl Genomes pan-taxonomic homologs.
##    subprocess.call(['perl', ensemblGenomesHomologScript, ensemblGeneIDs, ensemblGenomesHomologData, ensemblGenomesTaxonomyMap])

    # Parse the information from the Ensembl database.
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
    
    # Parse the information from the Ensembl Genome database.
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
    
    # Output the parse homology data. Preference is given to data extracted from Ensembl Genome as it is (potentially)
    # more up to date.
    homologs = []
    ensemblGenomeHomologs = set(ensemblGenomeHomologDict.keys())
    for i in ensemblGenomeHomologs:
        homologs.append(ensemblGenomeHomologDict[i])
    for i in ensemblHomologDict.keys():
        if i in ensemblGenomeHomologs:
            continue
        else:
            homologs.append(ensemblHomologDict[i])
    
    utilities.list2file.main(homologs, ensemblParsedHomology)