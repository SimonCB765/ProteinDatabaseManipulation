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
    sequenceOfAllProteinsDict = dict([(i, uniprotDict[i][-1]) for i in uniprotDict])  # A dictionary indexed by UniProt accession with the value for each index being the sequence of the protein.

    #===========================================================================
    # Add the data to the database.
    #===========================================================================
    rowsToAdd = [uniprotDict[i] for i in uniprotDict]
    values = '(' + ('%s,' * len(uniprotData[0]))
    values = values[:-1] + ')'
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    cursor = mysql.tableINSERT(cursor, tableProteinInfo, values, rowsToAdd)
    mysql.closeConnection(conn, cursor)
    print '\tEntries added to the UniProt table: ', len(toAdd)

    #===========================================================================
    # Annotate the proteins which have just been added to the database.
    #===========================================================================
    # Calculate the number of low complexity regions.
    proteinFasta = folderSEG + '/TempSEGFasta.fasta'
    SEGOutput = folderSEG + '/TempSEGOutput.txt'
    calculate_low_complexity(sequenceOfAllProteinsDict, SEGExe, proteinFasta, SEGOutput, schemaProteins,
                             tableProteinInfo, databasePassword)
    os.remove(proteinFasta)
    os.remove(SEGOutput)

    # Calculate the number of pest motifs.
    proteinFasta = folderEpestfind + '/TempEpestfindFasta.fasta'
    epestfindOutput = folderEpestfind + '/TempEpestfindOutput.txt'
    calculate_pest_motif(sequenceOfAllProteinsDict, epestfindExe, proteinFasta, epestfindOutput, schemaProteins,
                         tableProteinInfo, databasePassword)
    os.remove(proteinFasta)
    os.remove(epestfindOutput)

    # Calculate simple sequence statistics.
    proteinFasta = folderPepstats + '/TempPepstatsFasta.fasta'
    pepstatsOutput = folderPepstats + '/TempPepstatsOutput.txt'
    calculate_sequence_stats(sequenceOfAllProteinsDict, pepstatsExe, proteinFasta, pepstatsOutput, schemaProteins,
                             tableProteinInfo, databasePassword)
    os.remove(proteinFasta)
    os.remove(pepstatsOutput)

    #===========================================================================================
    # Use BLAST to determine the pairwise sequence identity of all the proteins in the table.
    #===========================================================================================
    BLASTOutput = folderBLAST + '/ProcessedBLAST.txt'
    if os.path.isfile(BLASTOutput):
        os.remove(BLASTOutput)

    # Generate the two fasta file to make the BLAST database from.
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

    # BLAST the proteins.
    tempBlastDatabaseFolder = folderBLAST + '/TempAllProtDB'
    os.mkdir(tempBlastDatabaseFolder)
    tempBlastDatabase = tempBlastDatabaseFolder + '/AllProt'
    makeDBArgs = makeBLASTDatabaseExe + ' -in ' + allProteinsFasta + ' -out ' + tempBlastDatabase
    subprocess.call(makeDBArgs)
    allBLASTOutput = folderBLAST + '/AllBLASTOutput.txt'
    out = ' -out ' + allBLASTOutput
    db = ' -db ' + tempBlastDatabase
    for i in sequenceOfAllProteinsDict.keys():
        writeTo = open(tempQuery, 'w')
        writeTo.write('>' + i + '\n')
        writeTo.write(sequenceOfAllProteinsDict[i] + '\n')
        writeTo.close()
        query = ' -query ' + tempQuery
        argsPSI = (query + out + evalue + inclusionEThresh + numIterations + gapTrigger + numDescriptions +
                   numAlignments + dbsize + db + outputFormat + numThreads)
        subprocess.call(psiblastExe + argsPSI, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        parsers.parsePSIBLAST.main(allBLASTOutput, BLASTOutput)
    shutil.rmtree(tempBlastDatabaseFolder)
    os.remove(allBLASTOutput)
    os.remove(tempQuery)
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
        if not query == hit:
            tuplesToAdd.append(matchesFound[i]['Tuple'])

    if tuplesToAdd != []:
        values = '(' + ('%s,' * len(tuplesToAdd[0]))
        values = values[:-1] + ')'
        cursor.execute('TRUNCATE TABLE ' + tableBLASTResults)
        cursor = mysql.tableINSERT(cursor, tableBLASTResults, values, tuplesToAdd)
    mysql.closeConnection(conn, cursor)

def calculate_low_complexity(sequenceOfAllProteinsDict, SEGExe, proteinFasta, SEGOutput, schemaProteins,
                             tableProteinInfo, databasePassword):

    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    for i in sequenceOfAllProteinsDict.keys():
        # Run every protein in the table through segmasker.
        UPAcc = i
        seq = sequenceOfAllProteinsDict[i]

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

def calculate_pest_motif(sequenceOfAllProteinsDict, epestfindExe, proteinFasta, epestfindOutput, schemaProteins,
                         tableProteinInfo, databasePassword):

    # Connect to the specified schema.
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    for i in sequenceOfAllProteinsDict.keys():
        # Run every protein in the table through epestfind.
        UPAcc = i
        seq = sequenceOfAllProteinsDict[i]

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

def calculate_sequence_stats(sequenceOfAllProteinsDict, pepstatsExe, proteinFasta, pepstatsOutput, schemaProteins,
                             tableProteinInfo, databasePassword):

    # Create the lists of the different types of amino acids. trueAAs are the 20 amino acids that are coded for by the genetic code.
    aminoAcids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    trueAAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    numAA = len(aminoAcids)
    tinyAAs = ['A', 'C', 'G', 'S', 'T']
    smallAAs = ['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V']
    aliphaticAAs = ['I', 'L', 'V']
    aromaticAAs = ['F', 'H', 'W', 'Y']
    nonpolarAAs = ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']
    polarAAs = ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T']
    chargedAAs = ['D', 'E', 'H', 'K', 'R']
    basicAAs = ['H', 'K', 'R']
    negativelyCharged = ['D', 'E']
    positivelyCharged = ['H', 'K', 'R']
    # The hydrophobicity of different amino acid residues, as measured by the Kyte and Doolittle scale.
    hydro = {'A' : 1.8, 'C' : 2.5, 'D' : -3.5, 'E' : -3.5, 'F' : 2.8, 'G' : -0.4, 'H' : -3.2, 'I' : 4.5,
             'K' : -3.9, 'L' : 3.8, 'M' : 1.9, 'N' : -3.5, 'P' : -1.6, 'Q' : -3.5, 'R' : -4.5, 'S' : -0.8,
             'T' : -0.7, 'V' : 4.2, 'W' : -0.9, 'Y' : -1.3}

    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)

    #===========================================================================
    # Run pepstats on all the proteins.
    #===========================================================================
    pepstatsInput = open(proteinFasta, 'w')
    for i in sequenceOfAllProteinsDict.keys():
        UPAcc = i
        seq = sequenceOfAllProteinsDict[i]
        pepstatsInput.write('>' + UPAcc + '\n')
        pepstatsInput.write(seq + '\n')
    pepstatsInput.close()

    # Run Pepstats on the newly created fasta file.
    subprocess.call(pepstatsExe + ' -sequence ' + proteinFasta + ' -outfile ' + pepstatsOutput + ' -auto')

    # Parse the Pepstats output file to get the isoelectric point of the protein.
    pIDict = parsers.parsePepstats.main(pepstatsOutput)

    for i in sequenceOfAllProteinsDict.keys():
        UPAcc = i
        pI = pIDict[UPAcc]['pI']

        # Write the isoelectric point into the Isoelectric column of the protein being analysed.
        cursor = mysql.tableUPDATE(cursor, tableProteinInfo, 'Isoelectric=' + str(pI), 'UPAccession = \'' + UPAcc + '\'')

    #===========================================================================
    # Calculate the sequence statistics for the proteins.
    #===========================================================================
    for i in sequenceOfAllProteinsDict.keys():
        stats = [0.0]*numAA  # The summation of the number of each type of amino acid in the protein.
        UPAcc = i
        seq = sequenceOfAllProteinsDict[i]
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