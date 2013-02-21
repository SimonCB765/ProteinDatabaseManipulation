'''
Created on 27 Mar 2012

@author: Simonial
'''

def main(targetDict, nonTargetDict, outputDirectory):

    outputData = {'biological_process' : {'LevelOne' : {'Unlabelled' : {}, 'Positive' : {}}, 'LevelTwo' : {'Unlabelled' : {}, 'Positive' : {}}},
                  'cellular_component' : {'LevelOne' : {'Unlabelled' : {}, 'Positive' : {}}, 'LevelTwo' : {'Unlabelled' : {}, 'Positive' : {}}},
                  'molecular_function' : {'LevelOne' : {'Unlabelled' : {}, 'Positive' : {}}, 'LevelTwo' : {'Unlabelled' : {}, 'Positive' : {}}}
                  }

    # Generate the summary data for the Unlabelled observations in the dataset.
    mainIndex = 'Unlabelled'
    for i in nonTargetDict.keys():
        proteinGOTerms = nonTargetDict[i]
        for j in proteinGOTerms.keys():
            levelOne = proteinGOTerms[j]['LevelOne']
            for k in levelOne:
                if outputData[j]['LevelOne'][mainIndex].has_key(k):
                    outputData[j]['LevelOne'][mainIndex][k] += 1
                else:
                    outputData[j]['LevelOne'][mainIndex][k] = 1
            levelTwo = proteinGOTerms[j]['LevelTwo']
            for k in levelTwo:
                if outputData[j]['LevelTwo'][mainIndex].has_key(k):
                    outputData[j]['LevelTwo'][mainIndex][k] += 1
                else:
                    outputData[j]['LevelTwo'][mainIndex][k] = 1

    # Generate the summary data for the Positive observations in the dataset.
    mainIndex = 'Positive'
    for i in targetDict.keys():
        proteinGOTerms = targetDict[i]
        for j in proteinGOTerms.keys():
            levelOne = proteinGOTerms[j]['LevelOne']
            for k in levelOne:
                if outputData[j]['LevelOne'][mainIndex].has_key(k):
                    outputData[j]['LevelOne'][mainIndex][k] += 1
                else:
                    outputData[j]['LevelOne'][mainIndex][k] = 1
            levelTwo = proteinGOTerms[j]['LevelTwo']
            for k in levelTwo:
                if outputData[j]['LevelTwo'][mainIndex].has_key(k):
                    outputData[j]['LevelTwo'][mainIndex][k] += 1
                else:
                    outputData[j]['LevelTwo'][mainIndex][k] = 1

    for i in outputData.keys():
        for j in outputData[i].keys():
            allKeys = sorted(set(outputData[i][j]['Unlabelled'].keys()).union(outputData[i][j]['Positive'].keys()))
            normalisedOutputData = {'Unlabelled' : {}, 'Positive' : {}}
            for k in allKeys:
                if outputData[i][j]['Unlabelled'].has_key(k):
                    normalisedOutputData['Unlabelled'][k] = outputData[i][j]['Unlabelled'][k]
                else:
                    normalisedOutputData['Unlabelled'][k] = 0
                if outputData[i][j]['Positive'].has_key(k):
                    normalisedOutputData['Positive'][k] = outputData[i][j]['Positive'][k]
                else:
                    normalisedOutputData['Positive'][k] = 0
            outputLocation = outputDirectory + '/' + i.upper() + '_' + j + '.txt'
            writeOut = open(outputLocation, 'w')
            writeOut.write('Class\t' + '\t'.join([k.replace(' ', '_') for k in allKeys]) + '\n')
            writeOut.write('1\t' + '\t'.join([str(normalisedOutputData['Unlabelled'][k]) for k in allKeys]) + '\n')
            writeOut.write('2\t' + '\t'.join([str(normalisedOutputData['Positive'][k]) for k in allKeys]) + '\n')
            writeOut.close()