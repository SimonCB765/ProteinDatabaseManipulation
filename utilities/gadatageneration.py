'''
Created on 27 Feb 2012

@author: Simonial
'''

import os
import shutil
import random

def pulearning(sqlPositiveQueryResults, sqlUnlabelledQueryResults, columnHeadings, outputLocation, categoricalMappingDataLocation, columnDataLocation,
            ECDataLocation, subcellLocation, healthStateLocation, bodySiteLocation, developmentalStageLocation, missingValueCode='-999'):
    """Generates the dataset in the format required by the Java random forest package.

    sqlPositiveQueryResults @type - list
    sqlPositiveQueryResults @use  - A list of data tuples. One tuple for each positive observation.
    sqlUnlabelledQueryResults @type - list
    sqlUnlabelledQueryResults @use  - A list of data tuples. One tuple for each unlabelled observation.
    columnHeadings @type - list
    columnHeadings @use  - A list of the names of the dependent variables in the dataset.
    outputLocation @type - string
    outputLocation @use  - The location where the output dataset will be stored.
    categoricalMappingDataLocation @type - string
    categoricalMappingDataLocation @use  - Categorical data is mapped to the range 1:numberOfCategories. This location will contain the mapping from category names to the 1:numberOfCategories range.
    columnDataLocation @type - string
    columnDataLocation @use  - The location where the names of the dependent variables are stored.
    ECDataLocation @type - string
    ECDataLocation @use  - The location where the EC number information should be stored.
    subcellLocation @type - string
    subcellLocation @use  - The location where the subcellular location information should be stored.
    healthStateLocation @type - string
    healthStateLocation @use  - The location where the health state expression information should be stored.
    bodySiteLocation @type - string
    bodySiteLocation @use  - The location where the body site expression information should be stored.
    developmentalStageLocation @type - string
    developmentalStageLocation @use  - The location where the developmental stage expression information should be stored.
    missingValueCode @type - string
    missingValueCode @use  - The value to give to variable observations that are missing.
    """

    columnIndices = dict([(columnHeadings[i], i) for i in range(len(columnHeadings))])

    # Determine the columns that will not be used when training the random forest.
    unneededColumns = ['UPAccession', 'ECNumber', 'SubcellularLocation', 'TopologicalDomain',
                       'PredictedSubcellularLocation', 'AlphaHelices', 'BetaStrands', 'PredictedBetaSheets', 'PredictedAlphaHelices']
    # Define the categorical/discrete columns.
    discreteColumns = []
    #discreteColumns = ['PESTMotif', 'LowComplexity', 'OGlycosylation', 'NGlycosylation', 'Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine',
    #                   'SignalPeptide', 'TransmembraneHelices', '3Untranslated', '5Untranslated', 'SynonymousCoding', 'Paralogs', 'BinaryPPI',
    #                   'HalfLife', 'InstabilityIndex']
    # Define the dependant variables that will be present in the training dataset, but are not to be used (do it this way so that they can be re-introduced easily).
    markedOutColumns = ['HalfLife', 'InstabilityIndex']
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
            healthStates = {}  # Records relative health state expression level information for the protein.
            bodySites = {}  # Records relative body site expression level information for the protein.
            developmentStages = {}  # Records relative developmental stage expression level information for the protein.
            for k in columnHeadings:
                if k == 'UPAccession':
                    # If the column heading is UPAccession then record the name of the protein, and whether it is positive or unlabelled.
                    protein = j[columnIndices[k]]
                    dataDict[protein] = {}
                    dataDict[protein]['Classification'] = 'Positive' if i == sqlPositiveQueryResults else 'Unlabelled'
                elif k == 'ECNumber':
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
                    # Get the Unigene health state expression information.
                    healthStates[k] = int(j[columnIndices[k]])
                elif k[:3] == 'BS_':
                    # Get the Unigene body site expression information.
                    bodySites[k] = int(j[columnIndices[k]])
                elif k[:3] == 'DS_':
                    # Get the Unigene developmental stage expression information.
                    developmentStages[k] = int(j[columnIndices[k]])
                elif k == 'InstabilityIndex':
                    # Determine whether the protein is predicted to be stable.
                    dataDict[protein][k] = '1' if j[columnIndices[k]] < 40 else '0'
                elif k == 'HalfLife':
                    # For a few of the proteins the N-terminus is not one of the twenty amino acids that ProtPAram works for.
                    # In these cases the value of the half life has been set to -1. When generating the data for the Orange data mining package,
                    # a half life of -1 should be recorded as a blank. Assumed to be correct procedure by looking at this page:
                    # http://orange.biolab.si/doc/reference/Orange.data.formats/
                    dataDict[protein][k] = str(j[columnIndices[k]]) if j[columnIndices[k]] != -1 else missingValueCode
                elif k in ['OGlycosylation', 'NGlycosylation', 'Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine',
                           'SignalPeptide', 'TransmembraneHelices', 'PredictedAlphaHelices', 'PredictedBetaSheets']:
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
            # Calculate the expression percentages for each of the three divisions.
            percentExpressionLevels[protein] = {'HealthState' : {}, 'BodySite' : {}, 'DevelopmentalStage' : {}}

            # Record relative data about the health state expression levels for the protein.
            totalHealthStateExpr = float(sum([healthStates[k] for k in healthStates.keys()]))
            healthStates = dict([(k, healthStates[k]/totalHealthStateExpr if healthStates[k] != 0 else 0) for k in healthStates.keys()])
            sortedExpressions = sorted([healthStates[k] for k in healthStates.keys()], reverse=True)
            if max(sortedExpressions) != 0:
                expressionPercentages = {0.1 : 0, 0.2 : 0, 0.3 : 0, 0.4 : 0, 0.5 : 0, 0.6 : 0, 0.7 : 0, 0.8 : 0, 0.9 : 0}
                for k in expressionPercentages.keys():
                    total = 0
                    index = 0
                    while total < k:
                        total += sortedExpressions[index]
                        index += 1
                    expressionPercentages[k] = index
                expressionPercentages[1.0] = sum([1 if k != 0 else 0 for k in sortedExpressions])
                percentExpressionLevels[protein]['HealthState'] = expressionPercentages
            else:
                percentExpressionLevels[protein]['HealthState'] = 'NA'

            # Record relative data about the body site expression levels for the protein.
            totalBodySiteExpr = float(sum([bodySites[k] for k in bodySites.keys()]))
            bodySites = dict([(k, bodySites[k]/totalBodySiteExpr if bodySites[k] != 0 else 0) for k in bodySites.keys()])
            for k in bodySites.keys():
                dataDict[protein][k] = str(bodySites[k])
            sortedExpressions = sorted([bodySites[k] for k in bodySites.keys()], reverse=True)
            if max(sortedExpressions) != 0:
                expressionPercentages = {0.1 : 0, 0.2 : 0, 0.3 : 0, 0.4 : 0, 0.5 : 0, 0.6 : 0, 0.7 : 0, 0.8 : 0, 0.9 : 0}
                for k in expressionPercentages.keys():
                    total = 0
                    index = 0
                    while total < k:
                        total += sortedExpressions[index]
                        index += 1
                    expressionPercentages[k] = index
                expressionPercentages[1.0] = sum([1 if k != 0 else 0 for k in sortedExpressions])
                percentExpressionLevels[protein]['BodySite'] = expressionPercentages
            else:
                percentExpressionLevels[protein]['BodySite'] = 'NA'

            # Record relative data about the developmental stage expression levels for the protein.
            totalDevelopmentStageExpr = float(sum([developmentStages[k] for k in developmentStages.keys()]))
            developmentStages = dict([(k, developmentStages[k]/totalDevelopmentStageExpr if developmentStages[k] != 0 else 0) for k in developmentStages.keys()])
            for k in developmentStages.keys():
                dataDict[protein][k] = str(developmentStages[k])
            sortedExpressions = sorted([developmentStages[k] for k in developmentStages.keys()], reverse=True)
            if max(sortedExpressions) != 0:
                expressionPercentages = {0.1 : 0, 0.2 : 0, 0.3 : 0, 0.4 : 0, 0.5 : 0, 0.6 : 0, 0.7 : 0, 0.8 : 0, 0.9 : 0}
                for k in expressionPercentages.keys():
                    total = 0
                    index = 0
                    while total < k:
                        total += sortedExpressions[index]
                        index += 1
                    expressionPercentages[k] = index
                expressionPercentages[1.0] = sum([1 if k != 0 else 0 for k in sortedExpressions])
                percentExpressionLevels[protein]['DevelopmentalStage'] = expressionPercentages
            else:
                percentExpressionLevels[protein]['DevelopmentalStage'] = 'NA'

    # Determine all the dependant variables that will be used when training the random forest.
    allColumns = [i for i in columnHeadings if i in continuousColumns or i in discreteColumns] + responseColumn

    # Determine the number of categories for each discrete column.
    numberCategoriesMapping = dict([(i, set([])) for i in discreteColumns])
    for i in discreteColumns:
        for j in dataDict:
            numberCategoriesMapping[i].add(dataDict[j][i])
    # Determine how the values of the discrete columns should be mapped to the range 1:number of categories.
    categoryMapping = dict([(i, {}) for i in discreteColumns])
    for i in discreteColumns:
        orderedCategories = sorted(numberCategoriesMapping[i])
        numberCategories = range(1, len(orderedCategories) + 1)
        mappingDict = dict([(orderedCategories[j - 1], str(j)) for j in numberCategories])
        categoryMapping[i] = mappingDict
    # Map the values of the discrete columns to the range 1:number of categories, and write out the mappings.
    writeOutMappings = open(categoricalMappingDataLocation, 'w')
    for i in discreteColumns:
        mappingOutput = [i]
        for j in categoryMapping[i]:
            mappingOutput.append(j + '->' + categoryMapping[i][j])
        writeOutMappings.write('\t'.join(mappingOutput) + '\n')
        for j in dataDict:
            dataDict[j][i] = categoryMapping[i][dataDict[j][i]]
    writeOutMappings.close()

    # Write out the dataset in the format expected by the Java random forest implementation.
    writeOut = open(outputLocation, 'w')
    # Create the header lines.
    headerLineOne = []
    headerLineTwo = []
    headerLineThree = []
    for i in allColumns:
        headerLineOne.append(i)
        if i in responseColumn:
            headerLineTwo.append('r')
            headerLineThree.append('')
        elif i in markedOutColumns:
            headerLineTwo.append('x')
            headerLineThree.append('')
        elif i in discreteColumns:
            headerLineTwo.append('c')
            headerLineThree.append(str(len(numberCategoriesMapping[i])))
        elif i in continuousColumns:
            headerLineTwo.append('n')
            headerLineThree.append('')
    # Write out the header lines.
    writeOut.write('\t'.join(headerLineOne) + '\n')
    writeOut.write('\t'.join(headerLineTwo) + '\n')
    writeOut.write('\t'.join(headerLineThree) + '\n')
    for i in dataDict:
        proteinInfo = []
        for j in allColumns:
            proteinInfo.append(dataDict[i][j])
        writeOut.write('\t'.join(proteinInfo))
        writeOut.write('\n')
    writeOut.close()

    # Write out the dependent variable names.
    writeOut = open(columnDataLocation, 'w')
    columnOut = '\n'.join(allColumns)
    writeOut.write(columnOut)
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

    # Write out the expression information.
    percentsUsed = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    writeOutHealthState = open(healthStateLocation, 'w')
    writeOutHealthState.write('\t'.join([str(i) for i in percentsUsed]) + '\tClassification\n')
    writeOutBodySite = open(bodySiteLocation, 'w')
    writeOutBodySite.write('\t'.join([str(i) for i in percentsUsed]) + '\tClassification\n')
    writeOutDevelopmentalStage = open(developmentalStageLocation, 'w')
    writeOutDevelopmentalStage.write('\t'.join([str(i) for i in percentsUsed]) + '\tClassification\n')
    for i in percentExpressionLevels.keys():
        classification = dataDict[i]['Classification']

        healthStateData = percentExpressionLevels[i]['HealthState']
        if healthStateData != 'NA':
            healthStateData = '\t'.join([str(healthStateData[j]) for j in percentsUsed])
            writeOutHealthState.write(healthStateData + '\t' + classification + '\n')

        bodySiteData = percentExpressionLevels[i]['BodySite']
        if bodySiteData != 'NA':
            bodySiteData = '\t'.join([str(bodySiteData[j]) for j in percentsUsed])
            writeOutBodySite.write(bodySiteData + '\t' + classification + '\n')

        developmentalStageData = percentExpressionLevels[i]['DevelopmentalStage']
        if developmentalStageData != 'NA':
            developmentalStageData = '\t'.join([str(developmentalStageData[j]) for j in percentsUsed])
            writeOutDevelopmentalStage.write(developmentalStageData + '\t' + classification + '\n')

    writeOutHealthState.close()
    writeOutBodySite.close()
    writeOutDevelopmentalStage.close()
