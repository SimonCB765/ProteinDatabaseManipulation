'''
Created on 27 Feb 2012

@author: Simonial
'''

import os
import shutil
import random

def fortran(sqlTargetQueryResults, sqlNonTargetQueryResults, columnHeadings, outputLocation, columnDataLocation,
            ECDataLocation, subcellLocation, healthStateLocation, bodySiteLocation, developmentalStageLocation, missingValueCode='-999'):
    """Takes as input the results of an sql select query, the names of the columns, the location to write out the
    result and the location to write out the column headings in the order they appear in the data table.
    """

    columnIndices = dict([(columnHeadings[i], i) for i in range(len(columnHeadings))])
    
    unneededColumns = ['UPAccession', 'ECNumber', 'SubcellularLocation', 'TopologicalDomain',
                       'PredictedSubcellularLocation', 'AlphaHelices', 'BetaStrands', 'PredictedBetaSheets', 'PredictedAlphaHelices']
    discreteColumns = ['HalfLife', 'InstabilityIndex', 'Classification']
    continuousColumns = [i for i in columnHeadings if i not in unneededColumns and i not in discreteColumns and i[:3] != 'HS_']
    
    dataDict = {}
    ECNumberDict = {'1' : dict([(str(i), 0) for i in range(1,8)]), '2' : dict([(str(i), 0) for i in range(1,8)])}  # Primary EC numbers go from 1 to 6. The number 7 is used to represent proteins with no EC number.
    subcellMapping = {'chlo' : 'Chloroplast', 'chlo_mito' : 'Chloroplast_Mitochondria', 'cysk' : 'Cytoskeleton',
                      'cysk_plas' : 'Cytoskeleton_Plasma membrane', 'cyto' : 'Cytosol', 'cyto_E.R.' : 'Cytosol_ER',
                      'cyto_mito' : 'Cytosol_Mitochondria', 'cyto_nucl' : 'Cytosol_Nulceus', 'cyto_pero' : 'Cytosol_Peroxisome',
                      'cyto_plas' : 'Cytosol_Plasma membrane', 'E.R.' : 'ER', 'E.R._golg' : 'ER', 'E.R._mito' : 'ER_Mitochondria',
                      'E.R._plas' : 'ER_Plasma membrane', 'E.R._vacu' : 'ER_Vacuolar membrane', 'extr' : 'Extracellular',
                      'extr_plas' : 'Extracellular_Plasma membrane', 'golg' : 'Golgi apparatus', 'golg_plas' : 'Golgi apparatus_Plasma membrane',
                      'lyso' : 'Lysosome', 'mito' : 'Mitochondria', 'mito_nucl' : 'Mitochondria_Nulceus', 'mito_pero' : 'Mitochondria_Peroxisome',
                      'mito_plas' : 'Mitochondria_Plasma membrane', 'nucl' : 'Nulceus', 'nucl_plas' : 'Nulceus_Plasma membrane',
                      'pero' : 'Peroxisome', 'plas' : 'Plasma membrane', 'vacu' : 'Vacuolar membrane', 'NoPrediction' : 'NoPrediction'}
    subcellLocDict = {'1' : dict([(subcellMapping[i], 0) for i in subcellMapping.keys()]),
                      '2' : dict([(subcellMapping[i], 0) for i in subcellMapping.keys()])
                      }
    percentExpressionLevels = {}
    for i in [sqlTargetQueryResults, sqlNonTargetQueryResults]:
        for j in i:
            healthStates = {}
            bodySites = {}
            developmentStages = {}
            for k in columnHeadings:
                if k == 'UPAccession':
                    # If the column heading is UPAccession then record the name of the protein, and whether it is
                    # a target.
                    protein = j[columnIndices[k]]
                    dataDict[protein] = {}
                    dataDict[protein]['Classification'] = '2' if i == sqlTargetQueryResults else '1'
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

            totalHealthStateExpr = float(sum([healthStates[k] for k in healthStates.keys()]))
            healthStates = dict([(k, healthStates[k]/totalHealthStateExpr if healthStates[k] != 0 else 0) for k in healthStates.keys()])
##            for k in healthStates.keys():
##                dataDict[protein][k] = str(healthStates[k])
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

    allColumns = continuousColumns + discreteColumns
    
    writeOut = open(outputLocation, 'w')
    for i in dataDict:
        proteinInfo = []
        for j in allColumns:
            proteinInfo.append(dataDict[i][j])
        proteinInfo = '\t'.join(proteinInfo)
        writeOut.write(proteinInfo)
        writeOut.write('\n')
    writeOut.close()

    writeOut = open(columnDataLocation, 'w')
    columnOut = '\n'.join(allColumns)
    writeOut.write(columnOut)
    writeOut.close()
    
    writeOut = open(ECDataLocation, 'w')
    ECNames = ECNumberDict[ECNumberDict.keys()[0]].keys()
    writeOut.write('Class\t' + '\t'.join(sorted(ECNames)) + '\n')
    print 'COMPARE THIS TO FILE:'
    for i in sorted(ECNumberDict.keys()):
        print i, ECNumberDict[i]
        outputValue = [i]
        for j in sorted(ECNames):
            outputValue.append(str(ECNumberDict[i][j]))
        writeOut.write('\t'.join(outputValue) + '\n')
    writeOut.close()
    
    writeOut = open(subcellLocation, 'w')
    subcellNames = [subcellMapping[i] for i in subcellMapping.keys()]
    writeOut.write('Class\t' + '\t'.join([i.replace(' ', '_') for i in subcellNames]) + '\n')
    print 'COMPARE THIS TO FILE:'
    for i in sorted(subcellLocDict.keys()):
        print i, subcellLocDict[i]
        outputValue = [i]
        for j in subcellNames:
            outputValue.append(str(subcellLocDict[i][j]))
        writeOut.write('\t'.join(outputValue) + '\n')
    writeOut.close()

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

def fortran_split(fortranStyleData, classes, splits, outputLocation):
    """
    assumes fortranStyleData has the classes as the last column
    the sampling is stratified
    """

    outputLocation = outputLocation + '/DataForFortran'
    if os.path.isdir(outputLocation):
        shutil.rmtree(outputLocation)
    os.mkdir(outputLocation)

    classes = [str(i) for i in classes]

    # Read in the input data.
    inputData = dict([(i, set([])) for i in classes])
    readIn = open(fortranStyleData, 'r')
    for line in readIn:
        chunks = line.split()
        classification = chunks[-1]
        inputData[classification].add(line)
    readIn.close()

    # Shuffle the input data. This is done to randomise the choice of data in the splits.
    inputData = [list(inputData[i]) for i in classes]
    for i in inputData:
        random.shuffle(i)

    numberOfInstances = [len(i) for i in inputData]  # The number of instances of each class.
    instancesInEachSplit = [i / splits for i in numberOfInstances]  # The floor (integer division) of the number of instances in each split.
    splitsForClass = [range(0, numberOfInstances[i], instancesInEachSplit[i]) for i in range(len(classes))]  # Splits for each class.
    tempSplitData = [[inputData[i][j:j+instancesInEachSplit[i]] for j in splitsForClass[i]] for i in range(len(classes))]
    splitData = []  # Split data, ordered in the same order as the classes input parameter. Takes the form [classes[0] splits, classes[1] splits,...]. Where classes[i] splits is a list of data examples that belong in the given split.
    for i in tempSplitData:
        if len(i) == splits:
            splitData.append(i)
        else:
            # This occurs when the number of splits does not evenly split the number of examples (e.g. 5 splits for 13 examples). In this case there
            # will be splits+1 elements in i. The final element will be a list of length < splits. The entries in this final small list will be
            # ditributed throughout the other elements of i. For example, with 5 splits and 13 examples, the 6th split will have 3 examples in it.
            # One of the 3 examples will go to split 1, another to split 2 and the final to split 3.
            tempSplits = i[:-1]
            finalSplit = i[-1]
            lenOfFinalSplit = len(finalSplit)
            for j in range(lenOfFinalSplit):
                tempSplits[j].append(finalSplit[j])
            splitData.append(tempSplits)

    # Generate the non-split specific parameter information.
    numVariables = len(splitData[0][0][0].split()) - 1  # Subtract 1 to discount the class.
    largerClassIndex = numberOfInstances.index(max(numberOfInstances))
    weightings = ['1.0' if i == largerClassIndex else str(round(numberOfInstances[largerClassIndex] / numberOfInstances[i])) for i in range(len(numberOfInstances))]
    
    # Generate a sub-directory for each split.
    for i in range(splits):
        # i determines the split to use as the testing split.
        subDir = outputLocation + '/' + str(i)
        try:
            os.mkdir(subDir)
        except:
            pass

        testFile = subDir + '/data.test'
        numTestingExamples = 0
        writeTest = open(testFile, 'w')
        for j in splitData:
            numTestingExamples += len(j[i])
            for k in j[i]:
                writeTest.write(k)
        writeTest.close()

        trainingFile = subDir + '/data.train'
        numTrainingExamples = 0
        writeTraining = open(trainingFile, 'w')
        for j in splitData:
            tempTrainingData = j[:i] + j[i+1:]
            trainingData = []
            for k in tempTrainingData:
                trainingData.extend(k)
            numTrainingExamples += len(trainingData)
            for k in trainingData:
                writeTraining.write(k)
        writeTraining.close()

        parameterFile = subDir + '/Parameters.txt'
        writeParam = open(parameterFile, 'w')
        writeParam.write(str(numTrainingExamples) + '\t' + str(numTestingExamples) + '\t' + str(numVariables) + '\t' + '\t'.join(weightings) + '\n')
        writeParam.close()

def orange(sqlTargetQueryResults, sqlNonTargetQueryResults, columnHeadings, outputLocation):
    """Takes as input the results of an sql select query, the names of the columns, the location to write out the result.
    """

    columnIndices = dict([(columnHeadings[i], i) for i in range(len(columnHeadings))])
    
    unneededColumns = ['UPAccession', 'ECNumber', 'SubcellularLocation', 'TopologicalDomain',
                       'PredictedSubcellularLocation', 'AlphaHelices', 'BetaStrands', '3Untranslated', '5Untranslated',
                       'NonSynonymousCoding', 'SynonymousCoding']
    discreteColumns = ['HalfLife', 'InstabilityIndex', 'Classification']
    continuousColumns = [i for i in columnHeadings if i not in unneededColumns and i not in discreteColumns]
    
    dataDict = {}
    for i in [sqlTargetQueryResults, sqlNonTargetQueryResults]:
        for j in i:
            for k in columnHeadings:
                if k == 'UPAccession':
                    # If the column heading is UPAccession then record the name of the protein, and whether it is
                    # a target.
                    protein = j[columnIndices[k]]
                    dataDict[protein] = {}
                    dataDict[protein]['Classification'] = 'Target' if i == sqlTargetQueryResults else 'NonTarget'
                elif k in unneededColumns:
                    # If the heading is one of these (except for UPAccession which is handeled separately), then it
                    # does not need to be recorded for the purpose of the genetic algorithm.
                    pass
                elif k == 'InstabilityIndex':
                    # Determine whether the protein is predicted to be stable.
                    dataDict[protein][k] = '1' if j[columnIndices[k]] < 40 else '0'
                elif k == 'HalfLife':
                    # For a few of the proteins the N-terminus is not one of the twenty amino acids that ProtPAram works for.
                    # In these cases the value of the half life has een set to -1. When generating the data for the Orange data mining package,
                    # a hlaf life of -1 should be recorded as a blank. Assumed to be correct procedure by looking at this page:
                    # http://orange.biolab.si/doc/reference/Orange.data.formats/
                    dataDict[protein][k] = str(j[columnIndices[k]]) if j[columnIndices[k]] != -1 else '?'
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

    attributeTypes = ['continuous'] * len(continuousColumns)
    attributeTypes.extend(['discrete'] * len(discreteColumns))
    classIdentification = [''] * (len(attributeTypes) - 1)
    classIdentification.append('class')  # The last column is the class information.
    allColumns = continuousColumns + discreteColumns
    headingsToWrite = '\t'.join(allColumns)
    attributesToWrite = '\t'.join(attributeTypes)
    classToWrite = '\t'.join(classIdentification)
    
    writeOut = open(outputLocation, 'w')
    writeOut.write(headingsToWrite)
    writeOut.write('\n')
    writeOut.write(attributesToWrite)
    writeOut.write('\n')
    writeOut.write(classToWrite)
    writeOut.write('\n')
    for i in dataDict:
        proteinInfo = []
        for j in allColumns:
            proteinInfo.append(dataDict[i][j])
        proteinInfo = '\t'.join(proteinInfo)
        writeOut.write(proteinInfo)
        writeOut.write('\n')
    writeOut.close()