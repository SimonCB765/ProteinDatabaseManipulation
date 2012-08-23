'''
Created on 3 Apr 2012

@author: Simonial
'''

import subprocess
import re

import utilities.MySQLaccess as mysql

def main(databasePassword, schemaPharmGKB, pharmGKBDiseases, pharmGKBDrugs, pharmGKBGenes, pharmGKBRelationships,
         pharmGKBInitialisation, MYSQLBIN):
    
    # Define the tables needed.
    tableDisease = schemaPharmGKB + '.disease'
    tableDrug = schemaPharmGKB + '.drug'
    tableGene = schemaPharmGKB + '.gene'
    tableDrug2Drug = schemaPharmGKB + '.drug2drug'
    tableDrug2DrugClass = schemaPharmGKB + '.drug2drugclass'
    tableDrug2Disease = schemaPharmGKB + '.drug2disease'
    tableDrug2Gene = schemaPharmGKB + '.drug2gene'
    tableDrugClass2Disease = schemaPharmGKB + '.drugclass2disease'
    tableDrugClass2Gene = schemaPharmGKB + '.drugclass2gene'
    tableGene2Disease = schemaPharmGKB + '.gene2disease'
    tableGene2Gene = schemaPharmGKB + '.gene2gene'
    tableDisease2MeSH = schemaPharmGKB + '.disease2mesh'
    tableDisease2SnoMED = schemaPharmGKB + '.disease2snomed'
    tableDisease2UMLS = schemaPharmGKB + '.disease2umls'
    
    # Create the schema and the tables.
    subprocess.call('mysql.exe -u root -p' + databasePassword + ' mysql < ' + pharmGKBInitialisation, shell=True, cwd=MYSQLBIN)
    
    # Parse the files that contain the data to enter into the database.    
    diseaseTuples = []
    meshTuples = []
    snomedTuples = []
    umlsTuples = []
    readIn = open(pharmGKBDiseases, 'r')
    readIn.readline()  # Strip off the header of the file.
    for line in readIn:
        chunks = line.split('\t')
        diseaseID = chunks[0]
        diseaseName = chunks[1]
        diseaseTuples.append(tuple([diseaseID, diseaseName]))
        externalVocab = chunks[4]
        mesh = re.findall('(?<=MeSH:)D[0-9]+', externalVocab)
        meshTuples.extend([tuple([diseaseID, i]) for i in mesh])
        snomed = re.findall('(?<=SnoMedCT:)[0-9]+', externalVocab)
        snomedTuples.extend([tuple([diseaseID, i]) for i in snomed])
        umls = re.findall('(?<=UMLS:)C[0-9]+', externalVocab)
        umlsTuples.extend([tuple([diseaseID, i]) for i in umls])
    readIn.close()
    
    drugTuples = []
    readIn = open(pharmGKBDrugs, 'r')
    readIn.readline()  # Strip off the header of the file.
    for line in readIn:
        chunks = line.split('\t')
        drugID = chunks[0]
        drugName = chunks[1]
        compoundType = chunks[5]
        if compoundType == 'Drug/Small Molecule':
            compoundType = 'Drug'
        xrefs = chunks[6]
        drugBank = re.findall('(?<=drugBank:)DB[0-9]+', xrefs)
        drugBank = drugBank[0] if drugBank != [] else 'NA'
        pubchemCID = re.findall('(?<=pubChemCompound:)[0-9]+', xrefs)
        pubchemCID = pubchemCID[0] if pubchemCID != [] else 'NA'
        drugTuples.append(tuple([drugID, drugName, compoundType, drugBank, pubchemCID]))
    readIn.close()
    
    geneTuples = []
    readIn = open(pharmGKBGenes, 'r')
    readIn.readline()  # Strip off the header of the file.
    for line in readIn:
        chunks = line.split('\t')
        geneID = chunks[0]
        entrezID = chunks[1]
        if entrezID == '':
            entrezID = '0'
        geneTuples.append(tuple([geneID, entrezID]))
    readIn.close()
    
    drug2drugTuples = []
    drug2drugClassTuples = {}
    drug2diseaseTuples = []
    drugClass2diseaseTuples = []
    drug2geneTuples = {}
    drugClass2geneTuples = {}
    gene2diseaseTuples = {}
    gene2geneTuples = {}
    readIn = open(pharmGKBRelationships, 'r')
    readIn.readline()  # Strip off the header of the file.
    for line in readIn:
        chunks = line.split('\t')
        entityOne = chunks[0].split(':')
        entityTwo = chunks[2].split(':')
        evidence = chunks[5]
        publication = 'Y' if 'Publication' in evidence else 'N'
        pathway = 'Y' if 'Pathway' in evidence else 'N'
        variant = 'Y' if 'Variant' in evidence else 'N'
        if entityOne[0] == 'Drug':
            if entityTwo[0] == 'Drug':
                drug2drugTuples.append(tuple([entityOne[1], entityTwo[1], pathway, publication, variant]))
            elif entityTwo[0] == 'Drug Class':
                key = tuple([entityOne[1], entityTwo[1]])
                if drug2drugClassTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drug2drugClassTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drug2drugClassTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drug2drugClassTuples[key][4] == 'Y' else 'N'
                    drug2drugClassTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drug2drugClassTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Gene':
                key = tuple([entityOne[1], entityTwo[1]])
                if drug2geneTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drug2geneTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drug2geneTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drug2geneTuples[key][4] == 'Y' else 'N'
                    drug2geneTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drug2geneTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Disease':
                drug2diseaseTuples.append(tuple([entityOne[1], entityTwo[1], pathway, publication, variant]))
        elif entityOne[0] == 'Gene':
            if entityTwo[0] == 'Drug':
                key = tuple([entityTwo[1], entityOne[1]])
                if drug2geneTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drug2geneTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drug2geneTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drug2geneTuples[key][4] == 'Y' else 'N'
                    drug2geneTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drug2geneTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Drug Class':
                key = tuple([entityTwo[1], entityOne[1]])
                if drugClass2geneTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drugClass2geneTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drugClass2geneTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drugClass2geneTuples[key][4] == 'Y' else 'N'
                    drugClass2geneTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drugClass2geneTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Disease':
                key = tuple([entityOne[1], entityTwo[1]])
                if gene2diseaseTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or gene2diseaseTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or gene2diseaseTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or gene2diseaseTuples[key][4] == 'Y' else 'N'
                    gene2diseaseTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    gene2diseaseTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Gene':
                key = tuple(sorted([entityOne[1], entityTwo[1]]))
                if gene2geneTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or gene2geneTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or gene2geneTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or gene2geneTuples[key][4] == 'Y' else 'N'
                    gene2geneTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    gene2geneTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
        elif entityOne[0] == 'Disease':
            if entityTwo[0] == 'Drug':
                drug2diseaseTuples.append(tuple([entityTwo[1], entityOne[1], pathway, publication, variant]))
            elif entityTwo[0] == 'Drug Class':
                drugClass2diseaseTuples.append(tuple([entityTwo[1], entityOne[1], pathway, publication, variant]))
            elif entityTwo[0] == 'Gene':
                key = tuple([entityTwo[1], entityOne[1]])
                if gene2diseaseTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or gene2diseaseTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or gene2diseaseTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or gene2diseaseTuples[key][4] == 'Y' else 'N'
                    gene2diseaseTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    gene2diseaseTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
        elif entityOne[0] == 'Drug Class':
            if entityTwo[0] == 'Drug':
                key = tuple([entityTwo[1], entityOne[1]])
                if drug2drugClassTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drug2drugClassTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drug2drugClassTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drug2drugClassTuples[key][4] == 'Y' else 'N'
                    drug2drugClassTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drug2drugClassTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
            elif entityTwo[0] == 'Disease':
                drugClass2diseaseTuples.append(tuple([entityOne[1], entityTwo[1], pathway, publication, variant]))
            elif entityTwo[0] == 'Gene':
                key = tuple([entityOne[1], entityTwo[1]])
                if drugClass2geneTuples.has_key(key):
                    pathEvidence = 'Y' if pathway == 'Y' or drugClass2geneTuples[key][2] == 'Y' else 'N'
                    pubEvidence = 'Y' if publication == 'Y' or drugClass2geneTuples[key][3] == 'Y' else 'N'
                    varEvidence = 'Y' if variant == 'Y' or drugClass2geneTuples[key][4] == 'Y' else 'N'
                    drugClass2geneTuples[key] = tuple([key[0], key[1], pathEvidence, pubEvidence, varEvidence])
                else:
                    drugClass2geneTuples[key] = tuple([key[0], key[1], pathway, publication, variant])
    readIn.close()
    
    # Enter the data into the database.
    conn, cursor = mysql.openConnection(databasePassword, schemaPharmGKB)
    
    values = '(' + ('%s,' * len(diseaseTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDisease, values, diseaseTuples)
    
    values = '(' + ('%s,' * len(meshTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDisease2MeSH, values, meshTuples)
    
    values = '(' + ('%s,' * len(snomedTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDisease2SnoMED, values, snomedTuples)
    
    values = '(' + ('%s,' * len(umlsTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDisease2UMLS, values, umlsTuples)    
    
    values = '(' + ('%s,' * len(drugTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrug, values, drugTuples)
    
    values = '(' + ('%s,' * len(geneTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableGene, values, geneTuples)

    values = '(' + ('%s,' * len(drug2drugTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrug2Drug, values, drug2drugTuples)

    drug2drugClassTuples = [drug2drugClassTuples[i] for i in drug2drugClassTuples.keys()]
    values = '(' + ('%s,' * len(drug2drugClassTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrug2DrugClass, values, drug2drugClassTuples)

    values = '(' + ('%s,' * len(drug2diseaseTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrug2Disease, values, drug2diseaseTuples)

    values = '(' + ('%s,' * len(drugClass2diseaseTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrugClass2Disease, values, drugClass2diseaseTuples)

    drug2geneTuples = [drug2geneTuples[i] for i in drug2geneTuples.keys()]
    values = '(' + ('%s,' * len(drug2geneTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrug2Gene, values, drug2geneTuples)

    drugClass2geneTuples = [drugClass2geneTuples[i] for i in drugClass2geneTuples.keys()]
    values = '(' + ('%s,' * len(drugClass2geneTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableDrugClass2Gene, values, drugClass2geneTuples)

    gene2diseaseTuples = [gene2diseaseTuples[i] for i in gene2diseaseTuples.keys()]
    values = '(' + ('%s,' * len(gene2diseaseTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableGene2Disease, values, gene2diseaseTuples)

    gene2geneTuples = [gene2geneTuples[i] for i in gene2geneTuples.keys()]
    values = '(' + ('%s,' * len(gene2geneTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableGene2Gene, values, gene2geneTuples)
    mysql.closeConnection(conn, cursor)