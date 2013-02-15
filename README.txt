When running on a Windows system it is best to start the script which will ultimately perform PSI-BLAST from the command prompt as this suppresses the constant popping up of command prompt windows everytime a new sequence is run.

Also should change the packet size allowed by MySQL (I changed it to 32M and stuck the my.ini file in C:\).

The output for the predictions must be called once before the input can occur. This is because the output function -l and -b make the predictions tables.


#################################################################
Parsers
#################################################################
For the output files, the name in brackets indicates the name of the variable in the controller.pythat corresponds to the file.
	FileName (variableName) - Info

parseBindingDB
parseCGC
	Takes the file that contains the cancer gene census data, and returns a file containing the processed CGC data.
	CGCParsedData (CGCParsed) - A tab separated (tsv) file, with three elements on each line.
		The first element is the NCBI gene ID of the gene.
		The second element is Y if the gene contains a cancer causing germline mutation, else it is N.
		The third element is Y if the gene contains a cancer causing somatic mutation, else it is N.
parseCGI
parseChEMBL
parseCOSMIC
parseDB
parseEnsembl
parseepestfind
parseGO
	Returns a file containing the parsed GO data.
	GOParsed (parsedGOOutput) - A file of 6-tuples, one on each line.
		The first element is the numerical identifier of the GO term.
		The second element is the name of the GO term.
		The third element is the category (biological_process, cellular component or molecular_function) that the term belongs to.
		The fourth element is all the paths from the term to the category it belongs to. The paths are separated from one another using semi-colons, and the elements of each path are separated from one another using '#'.
		The fifth element is all the level one terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the category in element three.
		The sixth element is all the level two terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the terms in element five.
parseHomologs
parsePepstats
parsePSIBLAST
parseSEG
parseTTD
	Takes the TTD target database data and the TTD drug data, and returns ont file of the UniProt accessions of the approved target proteins and one containing a mapping of approved targets to approved drugs.
	UPAccessions (TTDUPAccessions) - A file of UniProt accessions, one on each line.
		Each UP accession in the file is the target of an approved drug (as recorded by the TTD).
	TTDTarget2Drug (TTDTarget2Drug) - A tab separated (tsv) file, with three elements on each line.
		The first element is the TTD ID for the protein.
		The second element is a comma separated list of UniProt accessions that the first element is linked to.
		The third element is a semi-colon separated list of approved drugs that target the protein. For each drug a comma separated list of three elements is recorded.
			The first element is the name of the drug (as recorded by the TTD).
			The second element is the CAS number of the drug.
			The third element is the PubChem CID for the drug.
parseUG
	Takes the file that summaraiss the expression profile of ESTs in each cluster, and returns two processed files.
	UniGeneParsed (unigeneParsedOutput) - A file of tuples, one on each line.
		The first element of the tuple is the UniGene numerical Id of the cluster.
		The remaining elements of the tuple are records of the number of ESTs expressed for each expression option.
	UniGeneTotals (unigeneParsedTotals) - A file of 2-tuples, one on each line.
		The first element of the tuple is the body site, developmental stage or health site of the expression.
		The second element of the tuple is the number of ESTs for the given body site, developmental stage or health state over all cluster.
parseUP

#################################################################
Updaters
#################################################################
updateCancer
updateCOSMIC
updatedrug
updateEnsembl
updateGO
updatepathway
updatePPI
updatestability
updatetargetandredundancy
updateUG
updateUP
updatexref

#################################################################
Utilities
#################################################################
biomartquery
createfasta
file2list
gadatageneration
generateGOsummarydata
list2file
MySQLaccess
predictions
XMLparser