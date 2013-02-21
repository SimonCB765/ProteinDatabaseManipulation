'''
Created on 18 Oct 2011

@author: Simon Bull
'''

import subprocess
import os
import shutil

import parsers.parsePepstats
import parsers.parseSEG
import parsers.parseepestfind
import parsers.parsePSIBLAST

import utilities.file2list
import utilities.MySQLaccess as mysql

def main(UPProteinInfo, schemaProteins, tableProteinInfo, folderSEG, SEGExe, folderEpestfind, epestfindExe,
         folderPepstats, pepstatsExe, folderBLAST, psiblastExe, makeBLASTDatabaseExe, tableBLASTResults,
         databasePassword):

    #===========================================================================
    # Extract and format the parsed protein data.
    #===========================================================================
    uniprotData = utilities.file2list.main(UPProteinInfo)
    uniprotData = [eval(i) for i in uniprotData]
    uniprotDict = dict([(i[0], i) for i in uniprotData])

    #===========================================================================
    # Extract the protein information recorded in the database.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableSELECT(cursor, '*', tableProteinInfo)
    results = cursor.fetchall()
    mysql.closeConnection(conn, cursor)

    #===========================================================================
    # Compare the parsed data with the data recorded in the table.
    #===========================================================================
    sequenceIndex = 55
    columnIndices = [1, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 53, sequenceIndex]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor.execute('SHOW COLUMNS FROM ' + tableProteinInfo)
    columns = cursor.fetchall()
    mysql.closeConnection(conn, cursor)
    columns = [i[0] for i in columns]

    toRemove = []
    toUpdate = {}
    toAdd = uniprotDict.keys()
    for i in results:
        UPAcc = i[0]
        if uniprotDict.has_key(UPAcc):
            # If the key is in both the parsed file and the table, then it does not need to be added.
            toAdd.remove(i[0])
            # Compare the row from the table with the parsed file, to determine if the table needs updating.
            for j in columnIndices:
                if i[j] != uniprotDict[UPAcc][j]:
                    if not toUpdate.has_key(UPAcc):
                        toUpdate[UPAcc] = []
                    toUpdate[UPAcc].append(j)
        else:
            # If the key is in the table, but not in the parsed file, then the row needs to be removed.
            toRemove.append(i[0])
    values = '(' + ('%s,' * len(uniprotData[0]))
    values = values[:-1] + ')'

    # Record the proteins that have had their sequence changed/are newly added. These will need sequence properties
    # calculated, BLASTing and to have predictions made.
    sequenceChanged = [i for i in toUpdate.keys() if sequenceIndex in toUpdate[i]]
    sequenceChanged.extend([i for i in toAdd])
    sequenceChangedDict = dict([(i, uniprotDict[i][-1]) for i in sequenceChanged])
    sequenceNotChanged = [i[0] for i in results if i[0] not in sequenceChanged and i[0] not in toRemove]
    sequenceNotChangedDict = dict([(i, uniprotDict[i][-1]) for i in sequenceNotChanged])
    sequenceOfAllProteins = [i for i in set(sequenceChanged).union(sequenceNotChanged)]
    sequenceOfAllProteinsDict = dict([(i, uniprotDict[i][-1]) for i in sequenceOfAllProteins])

    #===========================================================================
    # Remove rows from the table that are not in the parsed file.
    #===========================================================================
    for i in toRemove:
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.rowDELETE(cursor, tableProteinInfo, 'UPAccession="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries removed from the UniProt table: ', len(toRemove)

    #===========================================================================
    # Update rows that have different values in the parsed file and the table.
    #===========================================================================
    for i in toUpdate.keys():
        toSet = []
        for j in toUpdate[i]:
            updateString = columns[j] + ' = "' + uniprotDict[i][j] + '"'
            toSet.append(updateString)
        toSet = ', '.join(toSet)
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, toSet, 'UPAccession="' + i + '"')
        mysql.closeConnection(conn, cursor)
    print '\tEntries updated in the UniProt table: ', len(toUpdate)

    #===========================================================================
    # Add rows which are not in the table, but are in the parsed file.
    #===========================================================================
    rowsToAdd = [uniprotDict[i] for i in toAdd]
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableINSERT(cursor, tableProteinInfo, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the UniProt table: ', len(toAdd)

    if len(sequenceChanged) > 0:
        # If some sequences have changed then the annotations and BLASTing can be performed.
        print '\tNow annotating added/altered proteins.'

        #===========================================================================
        # Annotate the proteins which have just been added to the table, or had
        # their sequence updated.
        #===========================================================================
        # Calculate the number of low complexity regions.
        proteinFasta = folderSEG + '/TempSEGFasta.fasta'
        SEGOutput = folderSEG + '/TempSEGOutput.txt'
        calculate_low_complexity(sequenceChangedDict, SEGExe, proteinFasta, SEGOutput, schemaProteins,
                                 tableProteinInfo, databasePassword)
        os.remove(proteinFasta)
        os.remove(SEGOutput)

        # Calculate the number of pest motifs.
        proteinFasta = folderEpestfind + '/TempEpestfindFasta.fasta'
        epestfindOutput = folderEpestfind + '/TempEpestfindOutput.txt'
        calculate_pest_motif(sequenceChangedDict, epestfindExe, proteinFasta, epestfindOutput, schemaProteins,
                             tableProteinInfo, databasePassword)
        os.remove(proteinFasta)
        os.remove(epestfindOutput)

        # Calculate simple sequence statistics.
        proteinFasta = folderPepstats + '/TempPepstatsFasta.fasta'
        pepstatsOutput = folderPepstats + '/TempPepstatsOutput.txt'
        calculate_sequence_stats(sequenceChangedDict, pepstatsExe, proteinFasta, pepstatsOutput, schemaProteins,
                                 tableProteinInfo, databasePassword)
        os.remove(proteinFasta)
        os.remove(pepstatsOutput)

        #===========================================================================
        # Use BLAST to determine the pairwise sequence identity of all the proteins
        # in the table. Not all sequence identities need to be calculated.
        #===========================================================================
        BLASTOutput = folderBLAST + '/ProcessedBLAST.txt'
        if os.path.isfile(BLASTOutput):
            os.remove(BLASTOutput)

        # Generate the two fasta files to make BLAST databases from.
        changedFasta = folderBLAST + '/TempChangeFasta.fasta'
        writeTo = open(changedFasta, 'w')
        for i in sequenceChangedDict.keys():
            writeTo.write('>' + i + '\n')
            writeTo.write(sequenceChangedDict[i] + '\n')
        writeTo.close()
        allProteinsFasta = folderBLAST + '/TempAllProteinsFasta.fasta'
        writeTo = open(allProteinsFasta, 'w')
        for i in sequenceOfAllProteinsDict.keys():
            writeTo.write('>' + i + '\n')
            writeTo.write(sequenceOfAllProteinsDict[i] + '\n')
        writeTo.close()

        # Set the non-unique PSI-BLAST parameters
        tempQuery = folderBLAST + '/TempQuery.fasta'
        evalue = ' -evalue 1'
        inclusionEThresh = ' -inclusion_ethresh 0.0001'
        numIterations = ' -num_iterations 3'
        gapTrigger = ' -gap_trigger 18'
        numDescriptions = ' -num_descriptions 10000'
        numAlignments = ' -num_alignments 10000'
        dbsize = ' -dbsize 0'
        outputFormat = ' -outfmt "7 qseqid sseqid pident length evalue"'
        numThreads = ' -num_threads 2'

        # Blast the new and updated proteins against all the proteins.
        tempBlastDatabaseFolder = folderBLAST + '/TempAllProtDB'
        os.mkdir(tempBlastDatabaseFolder)
        tempBlastDatabase = tempBlastDatabaseFolder + '/AllProt'
        makeDBArgs = makeBLASTDatabaseExe + ' -in ' + allProteinsFasta + ' -out ' + tempBlastDatabase
        subprocess.call(makeDBArgs)
        changedAgainstAllOutput = folderBLAST + '/ChangedAgainstAll.txt'
        out = ' -out ' + changedAgainstAllOutput
        db = ' -db ' + tempBlastDatabase
        for i in sequenceChangedDict.keys():
            writeTo = open(tempQuery, 'w')
            writeTo.write('>' + i + '\n')
            writeTo.write(sequenceChangedDict[i] + '\n')
            writeTo.close()
            query = ' -query ' + tempQuery
            argsPSI = (query + out + evalue + inclusionEThresh + numIterations + gapTrigger + numDescriptions +
                       numAlignments + dbsize + db + outputFormat + numThreads)
            subprocess.call(psiblastExe + argsPSI, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            parsers.parsePSIBLAST.main(changedAgainstAllOutput, BLASTOutput)
        shutil.rmtree(tempBlastDatabaseFolder)
        os.remove(changedAgainstAllOutput)
        os.remove(tempQuery)

        # Blast the non-updated proteins against the new and updated proteins, provided there are some unchanged ones.
        if len(sequenceNotChanged) > 0:
            tempBlastDatabaseFolder = folderBLAST + '/TempAddedAndChangedProtDB'
            os.mkdir(tempBlastDatabaseFolder)
            tempBlastDatabase = tempBlastDatabaseFolder + '/ChangedProt'
            makeDBArgs = makeBLASTDatabaseExe + ' -in ' + changedFasta + ' -out ' + tempBlastDatabase
            subprocess.call(makeDBArgs)
            sameAgainstChangedOutput = folderBLAST + '/SameAgainstChanged.txt'
            out = ' -out ' + sameAgainstChangedOutput
            db = ' -db ' + tempBlastDatabase
            for i in sequenceNotChangedDict.keys():
                writeTo = open(tempQuery, 'w')
                writeTo.write('>' + i + '\n')
                writeTo.write(sequenceNotChangedDict[i] + '\n')
                writeTo.close()
                query = ' -query ' + tempQuery
                argsPSI = (query + out + evalue + inclusionEThresh + numIterations + gapTrigger + numDescriptions +
                           numAlignments + dbsize + db + outputFormat + numThreads)
                subprocess.call(psiblastExe + argsPSI)
                parsers.parsePSIBLAST.main(sameAgainstChangedOutput, BLASTOutput)
            shutil.rmtree(tempBlastDatabaseFolder)
            os.remove(sameAgainstChangedOutput)
            os.remove(tempQuery)

        os.remove(changedFasta)
        os.remove(allProteinsFasta)

        # Enter the BLAST information into the blast results table.
        conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
        matchesFound = {}
        readIn = open(BLASTOutput, 'r')
        for line in readIn:
            chunks = line.split('\t')
            index = tuple(sorted([chunks[0], chunks[1]]))
            query = index[0]
            hit = index[1]
            similarity = float(chunks[2])
            length = int(chunks[3])
            eValue = float(chunks[4])
            if not matchesFound.has_key(index):
                # If this is the first time the (query, hit) pair has been found.
                matchesFound[index] = {}
                matchesFound[index]['Similarity'] = similarity
                matchesFound[index]['Tuple'] = tuple([query, hit, similarity, length, eValue])
            else:
                # If the (query, hit) pair has been found previously.
                if similarity < matchesFound[index]['Similarity']:
                    # If the similarity of the new (query, hit) pair is less than the previously found occurence of the
                    # pair, then discard this pair. This is because we are looking for the situation where the
                    # similarity is greatest, in order to know which pairs are redundant. Continue to next loop
                    # iteration as we do not want to update the database with this pair, or add this pair.
                    continue
                else:
                    matchesFound[index]['Similarity'] = similarity
                    matchesFound[index]['Tuple'] = tuple([query, hit, similarity, length, eValue])
        readIn.close()

        tuplesToAdd = []
        for i in matchesFound:
            query = matchesFound[i]['Tuple'][0]
            hit = matchesFound[i]['Tuple'][1]
            if query in toAdd or hit in toAdd:
                # If either of the query or hit protein is being added for the first time, then the tuple needs
                # inserting.
                tuplesToAdd.append(matchesFound[i]['Tuple'])
            else:
                # If neither protein is newly added, then either one of the proteins is already in the table and the
                # other protein is being updated, or both proteins are being updated. In this case the tuple in the
                # table needs to be updated.
                toSet = ('Similarity="' + str(matchesFound[i]['Tuple'][2]) + '", Length="' +
                         str(matchesFound[i]['Tuple'][3]) + '", EValue="' + str(matchesFound[i]['Tuple'][4]) + '"')
                cursor = mysql.tableUPDATE(cursor, tableBLASTResults, toSet,
                                           'ProteinA="' + query + '" AND ProteinB="' + hit + '"')
        if tuplesToAdd != []:
            values = '(' + ('%s,' * len(tuplesToAdd[0]))
            values = values[:-1] + ')'
            cursor = mysql.tableINSERT(cursor, tableBLASTResults, values, tuplesToAdd)
        mysql.closeConnection(conn, cursor)

def calculate_low_complexity(sequenceChangedDict, SEGExe, proteinFasta, SEGOutput, schemaProteins,
                             tableProteinInfo, databasePassword):

    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    for i in sequenceChangedDict.keys():
        # Run every protein in the table through segmasker.
        UPAcc = i
        seq = sequenceChangedDict[i]

        # Create a FASTA format file for the protein. This is the input format used for segmasker.
        SEGInput = open(proteinFasta, 'w')
        SEGInput.write('>' + UPAcc + '\n')
        SEGInput.write(seq)
        SEGInput.close()

        # Run segmasker on the fasta file just created.
        subprocess.call(SEGExe + ' -in ' + proteinFasta + ' -out ' + SEGOutput)

        # Parse the segmasker output file to determine if there are any low complexity regions.
        numLowComplexity = parsers.parseSEG.main(SEGOutput)

        # Write the number of low complexity regions into the LowComplexity column of the protein being analysed.
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'LowComplexity=' + str(numLowComplexity), 'UPAccession = \'' + UPAcc + '\'')


    mysql.closeConnection(conn, cursor)

def calculate_pest_motif(sequenceChangedDict, epestfindExe, proteinFasta, epestfindOutput, schemaProteins,
                         tableProteinInfo, databasePassword):

    # Connect to the specified schema.
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    for i in sequenceChangedDict.keys():
        # Run every protein in the table through epestfind.
        UPAcc = i
        seq = sequenceChangedDict[i]

        # Create a FASTA format file for the protein. This is the input format used for epestfind.
        epestfindInput = open(proteinFasta, 'w')
        epestfindInput.write('>' + UPAcc + '\n')
        epestfindInput.write(seq)
        epestfindInput.close()

        # Run epestfind on the fasta file just created.
        subprocess.call(epestfindExe + ' -sequence ' + proteinFasta + ' -outfile ' + epestfindOutput + ' -auto -window 10 -order score -graph none')

        # Parse the epestfind output file to determine if there is a valid PEST motif.
        motifPresent = parsers.parseepestfind.main(epestfindOutput)

        # Write the number of valid PEST motifs into the PESTMotif column of the protein being analysed.
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'PESTMotif=' + str(motifPresent), 'UPAccession = \'' + UPAcc + '\'')


    mysql.closeConnection(conn, cursor)

def calculate_sequence_stats(sequenceChangedDict, pepstatsExe, proteinFasta, pepstatsOutput, schemaProteins,
                             tableProteinInfo, databasePassword):

    # Create the lists of the different types of amino acids. trueAAs are the 20 amino acids that are coded for by the genetic code.
    aminoAcids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    trueAAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    numAA = len(aminoAcids)
    tinyAAs = ['A', 'C', 'G', 'S', 'T']
    smallAAs = ['A', 'B', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V']
    aliphaticAAs = ['I', 'L', 'V']
    aromaticAAs = ['F', 'H', 'W', 'Y']
    nonpolarAAs = ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']
    polarAAs = ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z']
    chargedAAs = ['B', 'D', 'E', 'H', 'K', 'R', 'Z']
    basicAAs = ['H', 'K', 'R']
    negativelyCharged = ['B', 'D', 'E', 'Z']
    positivelyCharged = ['H', 'K', 'O', 'R']
    # The hydrophobicity of different amino acid residues, as measured by the Kyte and Doolittle scale.
    hydro = {'A' : 1.8, 'C' : 2.5, 'D' : -3.5, 'E' : -3.5, 'F' : 2.8, 'G' : -0.4, 'H' : -3.2, 'I' : 4.5,
             'K' : -3.9, 'L' : 3.8, 'M' : 1.9, 'N' : -3.5, 'P' : -1.6, 'Q' : -3.5, 'R' : -4.5, 'S' : -0.8,
             'T' : -0.7, 'V' : 4.2, 'W' : -0.9, 'Y' : -1.3}

    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    #===========================================================================
    # Run pepstats on all the proteins.
    #===========================================================================
    pepstatsInput = open(proteinFasta, 'w')
    for i in sequenceChangedDict.keys():
        UPAcc = i
        seq = sequenceChangedDict[i]
        pepstatsInput.write('>' + UPAcc + '\n')
        pepstatsInput.write(seq + '\n')
    pepstatsInput.close()

    # Run Pepstats on the newly created fasta file.
    subprocess.call(pepstatsExe + ' -sequence ' + proteinFasta + ' -outfile ' + pepstatsOutput + ' -auto')

    # Parse the Pepstats output file to get the isoelectric point of the protein.
    pIDict = parsers.parsePepstats.main(pepstatsOutput)

    for i in sequenceChangedDict.keys():
        UPAcc = i
        pI = pIDict[UPAcc]['pI']

        # Write the isoelectric point into the Isoelectric column of the protein being analysed.
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Isoelectric=' + str(pI), 'UPAccession = \'' + UPAcc + '\'')

    #===========================================================================
    # Calculate the sequence statistics for the proteins.
    #===========================================================================
    for i in sequenceChangedDict.keys():
        stats = [0.0]*numAA  # The summation of the number of each type of amino acid in the protein.
        UPAcc = i
        seq = sequenceChangedDict[i]
        seqLen = len(seq)

        # Go through the amino acids in the sequence and sum up the different types.
        for aa in seq:
            index = aminoAcids.index(aa)
            stats[index] += 1

        # Compensate for the fact that not all amino acids recorded in the sequence will be from the 20 coded for by the genome.

##        # Remove O, U and X from the count of amino acids
##        O = stats[aminoAcids.index('O')]
##        U = stats[aminoAcids.index('U')]
##        X = stats[aminoAcids.index('X')]
##        seqLen = seqLen - O - U - X

        # B corresponds to asparagine (N) or aspartic acid (D)
        # Get the number of N and the number of D and treat a B as N/(N+D) asparagines and D/(N+D) aspartic acids
        B = stats[aminoAcids.index('B')]
        N = stats[aminoAcids.index('N')]
        D = stats[aminoAcids.index('D')]
        if B != 0:
            extraN = N / 2.#N / (N + D)
            extraD = D / 2.#D / (N + D)
            stats[aminoAcids.index('N')] += B * extraN
            stats[aminoAcids.index('D')] += B * extraD

        # J corresponds to leucine (L) or isoleucine (I)
        # Get the number of L and the number of I and treat a J as L/(L+I) leucines and I/(L+I) isoleucines
        J = stats[aminoAcids.index('J')]
        L = stats[aminoAcids.index('L')]
        I = stats[aminoAcids.index('I')]
        if J != 0:
            extraL = L / 2.#L / (L + I)
            extraI = I / 2.#I / (L + I)
            stats[aminoAcids.index('L')] += J * extraL
            stats[aminoAcids.index('I')] += J * extraI

        # Z corresponds to glutamine (Q) or glutamic acid (E)
        # Get the number of Q and the number of E and treat a Z as Q/(Q+E) glutamines and E/(Q+E) glutamic acids
        Z = stats[aminoAcids.index('Z')]
        Q = stats[aminoAcids.index('Q')]
        E = stats[aminoAcids.index('E')]
        if Z != 0:
            extraQ = Q / 2.#Q / (Q + E)
            extraE = E / 2.#E / (Q + E)
            stats[aminoAcids.index('Q')] += Z * extraQ
            stats[aminoAcids.index('E')] += Z * extraE

        hydroCalc = 0
        tinySum = 0
        smallSum = 0
        aliphaticSum = 0
        aromaticSum = 0
        nonpolarSum = 0
        polarSum = 0
        chargedSum = 0
        basicSum = 0
        negativelyChargedSum = 0
        positivelyChargedSum = 0
        for i in trueAAs:
            # For each of the 20 amino acids coded for by the genome, insert the amino acid frequency information into the table.
            insertValue = stats[aminoAcids.index(i)] / seqLen
            cursor = mysql.tableUPDATE(cursor, tableProteinInfo, i + '=' + str(insertValue), 'UPAccession = \'' + UPAcc + '\'')

            # Calculate the hydrophobicity information.
            hydroCalc = hydroCalc + (stats[aminoAcids.index(i)] * hydro[i])

            # Determine the number of tiny, small, etc etc amino acids.
            if i in tinyAAs:
                tinySum += stats[aminoAcids.index(i)]
            if i in smallAAs:
                smallSum += stats[aminoAcids.index(i)]
            if i in aliphaticAAs:
                aliphaticSum += stats[aminoAcids.index(i)]
            if i in aromaticAAs:
                aromaticSum += stats[aminoAcids.index(i)]
            if i in nonpolarAAs:
                nonpolarSum += stats[aminoAcids.index(i)]
            if i in polarAAs:
                polarSum += stats[aminoAcids.index(i)]
            if i in chargedAAs:
                chargedSum += stats[aminoAcids.index(i)]
            if i in basicAAs:
                basicSum += stats[aminoAcids.index(i)]
            if i in negativelyCharged:
                negativelyChargedSum += stats[aminoAcids.index(i)]
            if i in positivelyCharged:
                positivelyChargedSum += stats[aminoAcids.index(i)]

        # Calculate the mean hydrophobicity of the sequence.
        hydroCalc /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Hydrophobicity=' + str(hydroCalc), 'UPAccession = \'' + UPAcc + '\'')
        # Calculate the fraction of the sequence that is made up of each class of amino acids.
        tinySum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Tiny=' + str(tinySum), 'UPAccession = \'' + UPAcc + '\'')
        smallSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Small=' + str(smallSum), 'UPAccession = \'' + UPAcc + '\'')
        aliphaticSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Aliphatic=' + str(aliphaticSum), 'UPAccession = \'' + UPAcc + '\'')
        aromaticSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Aromatic=' + str(aromaticSum), 'UPAccession = \'' + UPAcc + '\'')
        nonpolarSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'NonPolar=' + str(nonpolarSum), 'UPAccession = \'' + UPAcc + '\'')
        polarSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Polar=' + str(polarSum), 'UPAccession = \'' + UPAcc + '\'')
        chargedSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Charged=' + str(chargedSum), 'UPAccession = \'' + UPAcc + '\'')
        basicSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Basic=' + str(basicSum), 'UPAccession = \'' + UPAcc + '\'')
        negativelyChargedSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'NegativelyCharged=' + str(negativelyChargedSum), 'UPAccession = \'' + UPAcc + '\'')
        positivelyChargedSum /= seqLen
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'PositivelyCharged=' + str(positivelyChargedSum), 'UPAccession = \'' + UPAcc + '\'')

    mysql.closeConnection(conn, cursor)