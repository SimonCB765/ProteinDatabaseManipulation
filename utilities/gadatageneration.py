'''
Created on 27 Feb 2012

@author: Simonial
'''

import os
import shutil
import random

def pulearning(sqlPositiveQueryResults, sqlUnlabelledQueryResults, columnHeadings, outputLocation,
            ECDataLocation, subcellLocation, healthStateLocation, bodySiteLocation, developmentalStageLocation):
    """Generates the dataset in the format required by the Java random forest package.

    sqlPositiveQueryResults @type - list
    sqlPositiveQueryResults @use  - A list of data tuples. One tuple for each positive observation.
    sqlUnlabelledQueryResults @type - list
    sqlUnlabelledQueryResults @use  - A list of data tuples. One tuple for each unlabelled observation.
    columnHeadings @type - list
    columnHeadings @use  - A list of the names of the dependent variables in the dataset.
    outputLocation @type - string
    outputLocation @use  - The location where the output dataset will be stored.
    ECDataLocation @type - string
    ECDataLocation @use  - The location where the EC number information should be stored.
    subcellLocation @type - string
    subcellLocation @use  - The location where the subcellular location information should be stored.
    """

    columnIndices = dict([(columnHeadings[i], i) for i in range(len(columnHeadings))])

    # Determine the columns that will not be used when training the random forest.
    unneededColumns = ['ECNumber', 'SubcellularLocation', 'TopologicalDomain', 'PredictedSubcellularLocation', 'PredictedBetaSheets', 'PredictedAlphaHelices']
    # Define the categorical/discrete columns.
    discreteColumns = []
    # Define the name of the response column.
    responseColumn = ['Classification']
    # Define the continuous valued columns.
    continuousColumns = [i for i in columnHeadings if i not in unneededColumns and i not in discreteColumns and i not in responseColumn and i[:3] != 'HS_']

    # Dictionary containing the entire dataset.
    dataDict = {}
    # Primary EC numbers go from 1 to 6. The number 7 is used to represent proteins with no EC number.
    ECNumberDict = {'Positive' : dict([(str(i), 0) for i in range(1,8)]), 'Unlabelled' : dict([(str(i), 0) for i in range(1,8)])}
    # The mapping from WoLFPsort subcellular location predictions to the descriptive names of the predictions.
    subcellMapping = {'chlo' : 'Chloroplast', 'chlo_mito' : 'Chloroplast_Mitochondria', 'cysk' : 'Cytoskeleton',
                      'cysk_plas' : 'Cytoskeleton_Plasma membrane', 'cyto' : 'Cytosol', 'cyto_E.R.' : 'Cytosol_ER',
                      'cyto_mito' : 'Cytosol_Mitochondria', 'cyto_nucl' : 'Cytosol_Nulceus', 'cyto_pero' : 'Cytosol_Peroxisome',
                      'cyto_plas' : 'Cytosol_Plasma membrane', 'E.R.' : 'ER', 'E.R._golg' : 'ER', 'E.R._mito' : 'ER_Mitochondria',
                      'E.R._plas' : 'ER_Plasma membrane', 'E.R._vacu' : 'ER_Vacuolar membrane', 'extr' : 'Extracellular',
                      'extr_plas' : 'Extracellular_Plasma membrane', 'golg' : 'Golgi apparatus', 'golg_plas' : 'Golgi apparatus_Plasma membrane',
                      'lyso' : 'Lysosome', 'mito' : 'Mitochondria', 'mito_nucl' : 'Mitochondria_Nulceus', 'mito_pero' : 'Mitochondria_Peroxisome',
                      'mito_plas' : 'Mitochondria_Plasma membrane', 'nucl' : 'Nulceus', 'nucl_plas' : 'Nulceus_Plasma membrane',
                      'pero' : 'Peroxisome', 'plas' : 'Plasma membrane', 'vacu' : 'Vacuolar membrane', 'NoPrediction' : 'NoPrediction'}
    # Dictionary containing the subcellular location prediction information.
    subcellLocDict = {'Positive' : dict([(subcellMapping[i], 0) for i in subcellMapping.keys()]),
                      'Unlabelled' : dict([(subcellMapping[i], 0) for i in subcellMapping.keys()])
                      }
    # Dictionary containing the expression level information.
    percentExpressionLevels = {}

    # Go through the observations and process the stored data into the format needed for the Java random forest.
    for i in [sqlPositiveQueryResults, sqlUnlabelledQueryResults]:
        for j in i:
            for k in columnHeadings:
                if k == 'UPAccession':
                    # If the column heading is UPAccession then record the name of the protein, and whether it is positive or unlabelled.
                    protein = j[columnIndices[k]]
                    dataDict[protein] = {}
                    dataDict[protein]['Classification'] = 'Positive' if i == sqlPositiveQueryResults else 'Unlabelled'

                # Determine how to handle the current column.
                if k == 'ECNumber':
                    classifiedAs = dataDict[protein]['Classification']
                    ECValue = j[columnIndices[k]]
                    if ECValue == 'NA':
                        primaryNumber = '7'
                    else:
                        primaryNumber = (ECValue.split('.'))[0]
                    ECNumberDict[classifiedAs][primaryNumber] += 1
                elif k == 'PredictedSubcellularLocation':
                    classifiedAs = dataDict[protein]['Classification']
                    subcellValue = j[columnIndices[k]].split(',')[0]
                    subcellLocDict[classifiedAs][subcellMapping[subcellValue]] += 1
                elif k in unneededColumns:
                    # If the heading is one of these (except for UPAccession which is handeled separately), then it
                    # does not need to be recorded for the purpose of the genetic algorithm.
                    pass
                elif k[:3] == 'HS_':
                    # Ignore the Unigene health state expression information.
                    pass
                elif k == 'InstabilityIndex':
                    # Determine whether the protein is predicted to be stable.
                    dataDict[protein][k] = '1' if j[columnIndices[k]] < 40 else '0'
                elif k == 'HalfLife':
                    # For a few of the proteins the N-terminus is not one of the twenty amino acids that ProtPAram works for.
                    # In these cases the value of the half life has been set to -1.
                    dataDict[protein][k] = str(j[columnIndices[k]])
                elif k in ['OGlycosylation', 'NGlycosylation', 'Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine',
                           'SignalPeptide', 'TransmembraneHelices', 'PredictedAlphaHelices', 'PredictedBetaSheets', 'Turns',
                           'AlphaHelices', 'BetaStrands']:
                    # All of these columns are either NA, or are a set of entries split up by ';'.
                    dataDict[protein][k] = '0' if j[columnIndices[k]] == 'NA' else str(j[columnIndices[k]].count(';') + 1)
                elif k == 'Sequence':
                    # Record the length of the sequence.
                    dataDict[protein]['Sequence'] = str(len(j[columnIndices[k]]))
                else:
                    # If the column heading is not in any of the previous categories, then the entry in the table
                    # for the column can serve as the actual value to store in the data dictionary. This is because
                    # the data for these columns in the dictionary is a numeric value.
                    dataDict[protein][k] = str(j[columnIndices[k]])

    # Determine all the dependant variables that will be used when training the random forest.
    allColumns = [i for i in columnHeadings if i in continuousColumns or i in discreteColumns] + responseColumn

    # Write out the dataset in the format expected by the Java random forest implementation.
    writeOut = open(outputLocation, 'w')
    writeOut.write('\t'.join(allColumns) + '\n')
    for i in dataDict:
        proteinInfo = []
        for j in allColumns:
            proteinInfo.append(dataDict[i][j])
        writeOut.write('\t'.join(proteinInfo))
        writeOut.write('\n')
    writeOut.close()

    # Write out the EC number information.
    writeOut = open(ECDataLocation, 'w')
    ECNames = ECNumberDict[ECNumberDict.keys()[0]].keys()
    writeOut.write('Class\t' + '\t'.join(sorted(ECNames)) + '\n')
    for i in sorted(ECNumberDict.keys()):
        outputValue = [i]
        for j in sorted(ECNames):
            outputValue.append(str(ECNumberDict[i][j]))
        writeOut.write('\t'.join(outputValue) + '\n')
    writeOut.close()

    # Write out the subcellular location information.
    writeOut = open(subcellLocation, 'w')
    subcellNames = [subcellMapping[i] for i in subcellMapping.keys()]
    writeOut.write('Class\t' + '\t'.join([i.replace(' ', '_') for i in subcellNames]) + '\n')
    for i in sorted(subcellLocDict.keys()):
        outputValue = [i]
        for j in subcellNames:
            outputValue.append(str(subcellLocDict[i][j]))
        writeOut.write('\t'.join(outputValue) + '\n')
    writeOut.close()