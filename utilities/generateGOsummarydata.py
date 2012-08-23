'''
Created on 27 Mar 2012

@author: Simonial
'''

def main(targetDict, nonTargetDict, outputDirectory):

    outputData = {'biological_process' : {'LevelOne' : {'1' : {}, '2' : {}}, 'LevelTwo' : {'1' : {}, '2' : {}}},
                  'cellular_component' : {'LevelOne' : {'1' : {}, '2' : {}}, 'LevelTwo' : {'1' : {}, '2' : {}}},
                  'molecular_function' : {'LevelOne' : {'1' : {}, '2' : {}}, 'LevelTwo' : {'1' : {}, '2' : {}}}
                  }

    mainIndex = '1'
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

    mainIndex = '2'
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
            allKeys = sorted(set(outputData[i][j]['1'].keys()).union(outputData[i][j]['2'].keys()))
            normalisedOutputData = {'1' : {}, '2' : {}}
            for k in allKeys:
                if outputData[i][j]['1'].has_key(k):
                    normalisedOutputData['1'][k] = outputData[i][j]['1'][k]
                else:
                    normalisedOutputData['1'][k] = 0
                if outputData[i][j]['2'].has_key(k):
                    normalisedOutputData['2'][k] = outputData[i][j]['2'][k]
                else:
                    normalisedOutputData['2'][k] = 0
            outputLocation = outputDirectory + '/' + i.upper() + '_' + j + '.txt'
            writeOut = open(outputLocation, 'w')
            writeOut.write('Class\t' + '\t'.join([k.replace(' ', '_') for k in allKeys]) + '\n')
            writeOut.write('1\t' + '\t'.join([str(normalisedOutputData['1'][k]) for k in allKeys]) + '\n')
            writeOut.write('2\t' + '\t'.join([str(normalisedOutputData['2'][k]) for k in allKeys]) + '\n')
            writeOut.close()