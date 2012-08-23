'''
Created on 10 Oct 2011

@author: Simon Bull
'''

import utilities.list2file

def main(parseTree, geneRIFDict, parsedEntrez):
    cancers, cancerTypes = cancer_genes(parseTree)
    cancerGenes = find_disease_genes(cancers, cancerTypes, geneRIFDict)
    
    # Turn the disease dictionaries into a list of tuples.
    insertGenes = [str(tuple([i] + cancerGenes[i])) for i in geneRIFDict.keys()]
    
    utilities.list2file.main(insertGenes, parsedEntrez)

def cancer_genes(parseTree):
    # Generate cancer lists.
    ovarian  = parseTree['DOID:2394'].collect_children_and_synonyms([])
    ovarian = set(ovarian)
    melanoma  = parseTree['DOID:1909'].collect_children_and_synonyms([])
    melanoma = set(melanoma)
    breast  = parseTree['DOID:1612'].collect_children_and_synonyms([])
    breast = set(breast)
    lung  = parseTree['DOID:1324'].collect_children_and_synonyms([])
    lung = set(lung)
    pancreatic  = parseTree['DOID:1793'].collect_children_and_synonyms([])
    pancreatic = set(pancreatic)
    liver  = parseTree['DOID:3571'].collect_children_and_synonyms([])
    liver = set(liver)
    colon  = parseTree['DOID:219'].collect_children_and_synonyms([])
    colon = set(colon)
    prostate  = parseTree['DOID:10283'].collect_children_and_synonyms([])
    prostate = set(prostate)
    testicular  = parseTree['DOID:2998'].collect_children_and_synonyms([])
    testicular = set(testicular)
    oesophageal  = parseTree['DOID:5041'].collect_children_and_synonyms([])
    oesophageal = set(oesophageal)
    stomach  = parseTree['DOID:10534'].collect_children_and_synonyms([])
    stomach = set(stomach)
    heart  = parseTree['DOID:117'].collect_children_and_synonyms([])
    heart = set(heart)
    oralCavity  = parseTree['DOID:8618'].collect_children_and_synonyms([])
    oralCavity = set(oralCavity)
    leukaemia  = parseTree['DOID:1240'].collect_children_and_synonyms([])
    leukaemia = set(leukaemia)
    lymphoma  = parseTree['DOID:0060058'].collect_children_and_synonyms([])
    lymphoma = set(lymphoma)
    intestinal  = parseTree['DOID:10155'].collect_children_and_synonyms([])
    intestinal = set(intestinal)
    brain = parseTree['DOID:1319'].collect_children_and_synonyms([])
    brain = set(brain)
    cancer = parseTree['DOID:162'].collect_children_and_synonyms([])
    cancer = set(cancer)
    # Order the cancer terms in a list. The order is such that terms later in the list (i.e. with bigger indices)
    # can not be subsets of the terms earlier in the list. Cancer is a superset of all the terms, and intestinal
    # is a superset of colon.
    cancers = [list(ovarian), list(melanoma), list(breast), list(lung), list(pancreatic), list(liver), list(colon),
               list(prostate), list(testicular), list(oesophageal), list(stomach), list(heart), list(oralCavity),
               list(leukaemia), list(lymphoma), list(intestinal), list(brain), list(cancer)]
    types = ['Ovarian', 'Melanoma', 'Breast', 'Lung', 'Pancreatic', 'Liver', 'Colon', 'Prostate', 'Testicular',
             'Oesophageal', 'Stomach', 'Heart', 'Oral Cavity', 'Leukaemia', 'Lymphoma', 'Intestinal', 'Brain', 'Cancer']
    
    return cancers, types

def find_disease_genes(diseases, diseaseTypes, geneRIFDict):
    diseaseGenes = dict([(i, ['N' for i in diseases]) for i in geneRIFDict.keys()])

    for i in geneRIFDict.keys():
        for j in geneRIFDict[i]:
            for k in range(len(diseases)):
                for l in diseases[k]:
                    if l in j:
                        diseaseGenes[i][k] = 'Y'
                        break

    return diseaseGenes