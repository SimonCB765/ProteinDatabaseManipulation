'''
Created on 10 Oct 2011

@author: Simon Bull
'''

import utilities.list2file

def main(unigeneProfiles, parsedUnigene, unigeneParsedTotals):
    # If new expressions are added to UniGene, or some that are not in the human list get added, then you can
    # simply add a column to the UniGene database table and add the value to the expression list in the
    # appropriate location.
    expression = ['embryoid body', 'blastocyst', 'fetus', 'neonate', 'infant', 'juvenile', 'adult', 'adrenal tumor',
                  'bladder carcinoma','breast (mammary gland) tumor', 'cervical tumor', 'chondrosarcoma',
                  'colorectal tumor', 'esophageal tumor', 'gastrointestinal tumor','germ cell tumor','glioma',
                  'head and neck tumor', 'kidney tumor', 'leukemia', 'liver tumor', 'lung tumor', 'lymphoma',
                  'non-neoplasia','normal','ovarian tumor', 'pancreatic tumor',
                  'primitive neuroectodermal tumor of the CNS', 'prostate cancer', 'retinoblastoma',
                  'skin tumor','soft tissue/muscle tissue tumor', 'uterine tumor', 'adipose tissue',
                  'adrenal gland', 'ascites', 'bladder', 'blood', 'bone', 'bone marrow','brain', 'cervix',
                  'connective tissue', 'ear', 'embryonic tissue', 'esophagus', 'eye', 'heart', 'intestine',
                  'kidney', 'larynx', 'liver', 'lung','lymph', 'lymph node', 'mammary gland', 'mouth', 'muscle',
                  'nerve', 'ovary', 'pancreas', 'parathyroid', 'pharynx', 'pituitary gland', 'placenta','prostate',
                  'salivary gland', 'skin', 'spleen', 'stomach', 'testis', 'thymus', 'thyroid', 'tonsil', 'trachea',
                  'umbilical cord', 'uterus','vascular']
    
    exprTotalIDs = ['DS_Embryoid_Body', 'DS_Blastocyst', 'DS_Fetus', 'DS_Neonate', 'DS_Infant', 'DS_Juvenile',
                  'DS_Adult', 'HS_Adrenal_Tumor', 'HS_Bladder_Carcinoma', 'HS_Breast_Mammary_Gland_Tumor',
                  'HS_Cervical_Tumor', 'HS_Chondrosarcoma', 'HS_Colorectal_Tumor', 'HS_Esophageal_Tumor',
                  'HS_Gastrointestinal_Tumor', 'HS_Germ_Cell_Tumor', 'HS_Glioma', 'HS_Head_And_Neck_Tumor',
                  'HS_Kidney_Tumor', 'HS_Leukemia_Tumor', 'HS_Liver_Tumor', 'HS_Lung_Tumor', 'HS_Lymphoma',
                  'HS_Non_neoplasia', 'HS_Normal', 'HS_Ovarian_Tumor', 'HS_Pancreatic_Tumor',
                  'HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS', 'HS_Prostate_Cancer', 'HS_Retinoblastoma',
                  'HS_Skin_Tumor', 'HS_Soft_Tissue_Muscle_Tissue_Tumor', 'HS_Uterine_Tumor',
                  'BS_Adipose_Tissue', 'BS_Adrenal_Gland', 'BS_Ascites', 'BS_Bladder', 'BS_Blood', 'BS_Bone',
                  'BS_Bone_Marrow', 'BS_Brain', 'BS_Cervix', 'BS_Connective_Tissue', 'BS_Ear',
                  'BS_Embryonic_Tissue', 'BS_Esophagus', 'BS_Eye', 'BS_Heart', 'BS_Intestine', 'BS_Kidney',
                  'BS_Larynx', 'BS_Liver', 'BS_Lung', 'BS_Lymph', 'BS_Lymph_Node', 'BS_Mammary_Gland',
                  'BS_Mouth', 'BS_Muscle', 'BS_Nerve', 'BS_Ovary', 'BS_Pancreas', 'BS_Parathyroid',
                  'BS_Pharynx', 'BS_Pituitary_Gland', 'BS_Placenta', 'BS_Prostate', 'BS_Salivary_Gland',
                  'BS_Skin', 'BS_Spleen', 'BS_Stomach', 'BS_Testis', 'BS_Thymus', 'BS_Thyroid', 'BS_Tonsil',
                  'BS_Trachea', 'BS_Umbilical_Cord', 'BS_Uterus', 'BS_Vascular']

    exprTotalValues = {}
    defaultExpression = [0] * len(expression)
    ID = None
    expressionProfiles = {}  # Indexed by the UniGene cluster ID
    
    readIn = open(unigeneProfiles, 'r')
    for line in readIn:
        if line[0] == '>':
            chunks = line.split(None, 1)
            chunks = chunks[-1].split('|')
            if not expressionProfiles.has_key(chunks[0]):
                ID = chunks[0]
                expressionProfiles[ID] = list(defaultExpression)
            continue
        line = line.strip()
        chunks = line.split('\t')
        insertLoc = expression.index(chunks[0])
        chunks = chunks[1].split(' / ')
        expressionProfiles[ID][insertLoc] = chunks[0]
        exprTotalValues[exprTotalIDs[insertLoc]] = chunks[1]
    readIn.close()
    
    # Create the lists of tuples for insertion.
    tupleList = [str(tuple([i] + expressionProfiles[i])) for i in expressionProfiles.keys()]
    utilities.list2file.main(tupleList, parsedUnigene)
    
    tupleList = [str(tuple([i, exprTotalValues[i]])) for i in exprTotalValues.keys()]
    utilities.list2file.main(tupleList, unigeneParsedTotals)