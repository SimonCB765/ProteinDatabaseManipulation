'''
Created on 27 Oct 2011

@author: Simon Bull
'''

import sys
import os
import shutil

import utilities.list2file
import utilities.file2list
import utilities.MySQLaccess as mysql

def main(tableToUse, outputDirectory, predictionDirection, seqsPerFile, maxSeqLength, schemaProteins,
         databasePassword, columnWithPredictions):

    if predictionDirection != 'IN':
        # Clear the output directory, and then remake it.
        if os.path.exists(outputDirectory):
            shutil.rmtree(outputDirectory)
        os.mkdir(outputDirectory)
    
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    if predictionDirection == 'OUTA':
        cursor = mysql.tableSELECT(cursor, 'UPAccession, Sequence', tableToUse)
        results = cursor.fetchall()
        results = ['>' + '\n'.join([j[0], j[1]]) for j in results if len(j[1]) <= maxSeqLength]
        length = len(results)
        numberOfSplits = length / seqsPerFile
        if length%seqsPerFile != 0:
            numberOfSplits += 1
        fileOutput = [ results[i * seqsPerFile : (i+1) * seqsPerFile] for i in range(numberOfSplits)]
        for i in range(len(fileOutput)):
            utilities.list2file.main(fileOutput[i], outputDirectory + '/' + str(i) + '.fasta')
    elif predictionDirection == 'OUTS':
        cursor = mysql.tableSELECT(cursor, 'UPAccession, Sequence', tableToUse, columnWithPredictions + '="NA"')
        results = cursor.fetchall()
        results = ['>' + '\n'.join([j[0], j[1]]) for j in results if len(j[1]) <= maxSeqLength]
        length = len(results)
        numberOfSplits = length / seqsPerFile
        if length%seqsPerFile != 0:
            numberOfSplits += 1
        fileOutput = [ results[i * seqsPerFile : (i+1) * seqsPerFile] for i in range(numberOfSplits)]
        for i in range(len(fileOutput)):
            utilities.list2file.main(fileOutput[i], outputDirectory + '/' + str(i) + '.fasta')
    elif predictionDirection == 'IN':
        predFolderFiles = os.listdir(outputDirectory)
        fastaFiles = [i for i in predFolderFiles if i.split('.')[1] == 'fasta']
        fastaFilesNoExtension = [i.split('.')[0] for i in fastaFiles]
        predFiles = [(i, i.split('.')[0] + '.fasta') for i in predFolderFiles if i.split('.')[0] in fastaFilesNoExtension and i not in fastaFiles]
        if len(predFiles) > len(fastaFiles):
            # There are too many prediction files.
            print 'There are too many prediction files, or too many files with the same name as the fasta files.'
            sys.exit()
        updatesToPerform = []
        for i in predFiles:
            predictionFile = outputDirectory + '/' + i[0]
            fastaFile = outputDirectory + '/' + i[1]
            # Gather the accessions of the proteins that were supposed to have been predicted.
            outputProteins = utilities.file2list.main(fastaFile)[::2]
            outputProteins = [i[1:] for i in outputProteins]
            # Gather the input data.
            inputProteins = utilities.file2list.main(predictionFile)
            updateTuples = wolf_prediction_entry(outputProteins, inputProteins)
            updatesToPerform.extend(updateTuples)
        updateQuery = 'UPDATE ' + tableToUse + '\n\tSET ' + columnWithPredictions + ' = CASE UPAccession\n'
        for i in updatesToPerform:
            updateQuery += '\t\tWHEN "' + i[0] + '" THEN "' + i[1] + '"\n'
        updateQuery += '\tEND\nWHERE UPAccession IN ("' + '","'.join([i[0] for i in updatesToPerform]) + '")'
        cursor.execute(updateQuery)
    mysql.closeConnection(conn, cursor)

def wolf_prediction_entry(outputProteinNames, inputPredictions):
    inputPredictions = [i for i in inputPredictions if i.split()[0] in outputProteinNames]
    inputPredictions = dict([(i.split(None, 2)[0], i.split(None, 2)[2]) for i in inputPredictions])
    inputProteins = inputPredictions.keys()

    tuplesToUpdate = []
    for i in outputProteinNames:
        if i in inputProteins:
            locations = inputPredictions[i].replace(' ', '')
            locations = locations.split(',')
            subcellLoc = []
            for j in locations:
                if j == '':
                    continue
                subloc = j.split(':')
                subcellLoc.append(subloc[0] + ',' + subloc[1])
            subcellLoc = ';'.join(subcellLoc)
            tuplesToUpdate.append(tuple([i, subcellLoc]))
        else:
            tuplesToUpdate.append(tuple([i, 'NoPrediction']))

    return tuplesToUpdate