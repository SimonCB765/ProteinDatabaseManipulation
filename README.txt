######################################

When running on a Windows system it is best to start the script which will ultimately perform PSI-BLAST from the command prompt as this suppresses the
constant popping up of command prompt windows everytime a new sequence is run.

Also should change the packet size allowed by MySQL (I changed it to 32M and stuck the my.ini file in C:\)

first the database must be remade then you can deal with the predictions
	additionally the output for the predictions must be called once before the input can occur
	this is because the output function -l and -b make the predictions tables

######################################

Databases:
======================================
Any entry that is not the UPID or the Sequence will default to NA. The UPID and Sequence must be present.

If you change what information is being written into the database then follow these steps:
1st: Change the information extracted from the UniProt human proteome (not essential but will speed up subsequent tasks)
2nd: Change the information recorded in the createDatabase.py script in the processUniProt procedure (essential)
3rd: Change the columns in the proteome database (tableNameProt) in the createDatabase.py script in the self procedure (essential)
4th: Change the columns in the target database (annotatedTargets) in the createDatabase.py script in the self procedure (essential)


Protein database columns are ordered as follows (top of the list corresponds to leftmost column while bottom corresponds to rightmost):
---------------------------------------
UniProt Accession number
Protein Name as recorded by UniProt
OMIM IDs linked to the protein
Next 20 columns are the different amino acids in alphabetic order (i.e. A is on the left and Y on the right)
	The value for each column is the fraction of the sequence which the amino acids makes up (i.e 0.055 in the A column means the sequence is 5.5% A)
Hydrophocity value
Isoelectric point
Fraction of amino acids in the sequence which are negatively charged
Fraction of amino acids in the sequence which are positively charged
Fraction of amino acids in the sequence which are basic
Fraction of amino acids in the sequence which are charged
Fraction of amino acids in the sequence which are polar
Fraction of amino acids in the sequence which are non-polar
Fraction of amino acids in the sequence which are aromatic
Fraction of amino acids in the sequence which are alophatic
Fraction of amino acids in the sequence which are small
Fraction of amino acids in the sequence which are tiny
Number of low complexity regions
number of PEST motifs
The type of protein or NA (type is one of G-protein coupled receptor, Ion channel, Kinase, Protease, Nuclear Receptor)
The EC number or NA if the number was not in UniProt
Glycosylation
	NX,in qualifier;OY,;CZ,in modifier; where N, O and C indicate the type of glycosylation and X, Y and Z are the positions in the sequence
	there is potentially an in qualifier between the , and the ; which means that there is an in qualifier (e.g. in isoform 2#)
Phosphorylation
	SER-X,;TYR-Y,in qualifier;THR-Z,; where SER, TYR and THR indicate the type of phosphorylation and X, Y and Z are the positions in the sequence
	there is potentially an in qualifier between the , and the ; which means that there is an in qualifier (e.g. in isoform 2#)
GO term IDs each separated by a ;
Subcellular location as recorded by UniProt (currently it is just free text)
Signal peptide cleavage
	Each entry takes the form A,B;C,D; where A and C are the start locations for two different cleavages and B and D the end locations
Transmembrane helices
	Each entry takes the form A,B;C,D; where A and C are the start locations for two different transmembrane helices and B and D the end locations
Turns as recorded by UniProt
Alpha helices as recorded by UniProt
Beta strands as recorded by UniProt
Isoforms as recorded by UniProt
	Each entry takes the form A,B;C,D; where A and C are the names by which UniProt refers to the isoforms
	B and D are the UniProt accession numbers for the isoforms
SNPs as recorded by UniProt
	Each entry takes the form aa1POSaa2,dbSNP,OMIM;
	dbSNP indicates the ID for the SNP as recorded in dbSNP, this may be blank if UniProt has not linked to the particular SNP
	OMIM is the OMIM ID for the particular SNP variant, this may be blank if there is no known OMIM ID linked to the variant
The sequence

Target database columns are ordered as follows (top of the list corresponds to leftmost column while bottom corresponds to rightmost):
---------------------------------------
All columns are the same except there are 4 new columns between the SNPs and the sequence. These columns record the databases where the targets were found.
If the target was found in the particular database then the entry will be a Y otherwise it will be an N.
From left to right the databases are:
ChEMBL (ChEMBL)
DrugBank (DB)
Therapeutic Target Database (TTD)
UniProt (UP)

goinfo database pathinfo table columns are ordered as follows (top of the list corresponds to leftmost column while bottom corresponds to rightmost):
---------------------------------------
GOTermIDs (in GO format == GO: followed by 7 digits)
GOName (name of the GO ID as a string)
GOType (which of the three superclasses the GO ID falls under, a single string)
GOPaths (the paths from the GO ID to the GOType, each term in a path is separated by a # and each path ends with a ;)
LevelOne (the different names of the GO terms one level below the superclass, each term ends with a ;)
LevelTwo (the different names of the GO terms two levels below the superclass, each term ends with a ;)

go database information
---------------------------------------
term2term stores is_a and part_of relationships in the form

id	relationship_type_id	term1	term2	complete

where the relationship is:	term2 is a/part of term1


graph_path stores all ancestor relationships in the form

id	term1	term2	relationship_type_id	distance	relation_distance

where the relationship is:	term1 is an ancestor of term2 with distance hops between them



GA
------------------------------------------
Set the number of attributes in the GA to one less than the total number (i.e. do NOT include the class attribute as this
	will be automatically added to every run)
	37 variables being used currently