def main(CGCData, CGCParsed):
    parsedCGC = parse_CGC(CGCData)

    writeTo = open(CGCParsed, 'w')
    for i in sorted(parsedCGC.keys()):
        writeTo.write(str(i) + '\t' + parsedCGC[i]['Germline'] + '\t' + parsedCGC[i]['Somatic'] + '\n')
    writeTo.close()


def parse_CGC(CGCData):
    readIn = open(CGCData, 'r')
    data = readIn.read()
    readIn.close()
    data = data.split('\r')[1:]  # For some reason the file is split on carriage returns, so the normal file reading method doesn't work. [1:] in order to skip the header line.
    parsedCGC = {}
    for line in data:
        chunks = line.split('\t')
        ncbiGene = int(chunks[2])
        somatic = 'Y' if 'yes' in chunks[5] else 'N'
        germline = 'Y' if 'yes' in chunks[6] else 'N'
        if parsedCGC.has_key(ncbiGene):
            somatic = 'Y' if somatic == 'Y' or parsedCGC[ncbiGene]['Somatic'] == 'Y' else 'N'
            germline = 'Y' if germline == 'Y' or parsedCGC[ncbiGene]['Germline'] == 'Y' else 'N'
        parsedCGC[ncbiGene] = {'Somatic' : somatic, 'Germline' : germline}

    return parsedCGC