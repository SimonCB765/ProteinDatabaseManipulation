'''
Created on 7 Oct 2011

@author: Simon Bull
'''

# Python distribution imports.
import sys
import getopt
import subprocess
import os
import shutil

# File/database dump loaders.
import loaders.chembl
import loaders.GO
import loaders.pharmgkb

# File/database parsers.
import parsers.parseGO
import parsers.parseTTD
import parsers.parseChEMBL
import parsers.parseDB
import parsers.parseUP
import parsers.parseCGC
import parsers.parseUG
import parsers.parseEnsembl
import parsers.parseHomologs
import parsers.parseBindingDB
import parsers.parseCOSMIC
import parsers.parseCGI

# Databse table updaters.
import updaters.updateCancer
import updaters.updateGO
import updaters.updateUG
import updaters.updateEnsembl
import updaters.updateUP
import updaters.updatetargetandredundancy
import updaters.updatexref
import updaters.updatePPI
import updaters.updatedrug
import updaters.updatepathway
import updaters.updatestability
import updaters.updateCOSMIC

# Utility scripts.
import utilities.biomartquery
import utilities.file2list
import utilities.list2file
import utilities.MySQLaccess as mysql
import utilities.predictions
import utilities.gadatageneration
import utilities.generateGOsummarydata

#===============================================================================
# Global Definitions
#===============================================================================
ROOT = 'E:/PhD/ProteinDatabase'
DATA = ROOT + '/Data'
SOURCE = ROOT + '/DatabaseManipulation'
DATABASEPASSWORD = 'root'
DATABASEUSERNAME = 'root'
MYSQLBIN = 'C:/Program Files/MySQL/MySQL Server 5.1/bin'

def print_help():
    print 'HELP'
    print 'See the README for further information.'

def main(args):

    #===========================================================================
    # Data File Locations.
    #===========================================================================
    # Folders containing the data.
    folderBindingDB = DATA + '/BindingDB'
    folderBiomart = folderEnsembl + '/Biomart'
    folderBLAST = DATA + '/BLAST'
    folderCancerTargets = DATA + '/CancerTargets'
    folderCGC = DATA + '/CancerGeneCensus'
    folderCGI = DATA + '/CancerGeneIndex'
    folderChEMBL = DATA + '/ChEMBL'
    folderCOSMIC = DATA + '/COSMIC'
    folderCulling = DATA + '/Culling'
    folderDB = DATA + '/DrugBank'
    folderEnsembl = DATA + '/Ensembl'
    folderEnsemblPerlAPI = folderEnsembl + '/PerlAPI'
    folderEpestfind = DATA + '/Epestfind'
    folderGAData = DATA + '/GeneticAlgorithmData'
    folderGO = DATA + '/GeneOntology'
    folderPathwayCommons = DATA + '/PathwayCommons'
    folderPepstats = DATA + '/Pepstats'
    folderPredictions = DATA + '/Predictions'
    folderSEG = DATA + '/SEG'
    folderTTD = DATA + '/TherapeuticTargetDatabase'
    folderUG = DATA + '/UniGene'
    folderUP = DATA + '/UniProt'

    # BindingDB files used.
    bindingSDF = folderBindingDB + '/BindingDB2D.sdf'
    binding2PubChem = folderBindingDB + '/BDB_cid.txt'
    bindingParsed = folderBindingDB + '/BindingDBParsed.txt'

    # Biomart files used.
    biomartScript = folderBiomart + '/biomartquery.pl'
    biomartQuery = folderBiomart + '/BiomartQuery.xml'

    # BLAST files used.
    psiblastExe = folderBLAST + '/psiblast.exe'
    makeBLASTDatabaseExe = folderBLAST + '/makeblastdb.exe'

    # File of all known FDA approved anti-cancer targets.
    cancerTargets = folderCancerTargets + '/CancerDrugs.txt'

    # Cancer Gene Census files used.
    CGCData = folderCGC + '/cancer_gene_census.tsv'  # A file containing the Cancer Gene Census data.
    CGCParsed = folderCGC + '/CGCParsedData.txt'  # A tab separated (tsv) file, with three elements on each line.
                                                  # The first element is the NCBI gene ID of the gene.
                                                  # The second element is Y if the gene contains a cancer causing germline mutation, else it is N.
                                                  # The third element is Y if the gene contains a cancer causing somatic mutation, else it is N.

    # Cancer Gene Index files used.
    CGIData = folderCGI + '/NCI_CancerIndex_allphases_disease.xml'
    CGIUPAccessions = folderCGI + '/CGIUPAccessions.txt'
    CGIHGNCIDs = folderCGI + '/CGIHGNCIDs.txt'

    # ChEMBL files used.
    completeChEMBLDatabase = folderChEMBL + '/ChEMBL.sql'  # The ChEMBL database.
    ChEMBLUPAccessions = folderChEMBL + '/UPAccessions.txt'
    ChEMBLCID = folderChEMBL + '/PubChemCIDs.txt'

    # COSMIC files used.
    cosmicData = folderCOSMIC + '/CosmicCompleteExport_v59_220512.tsv'
    cosmicParsedGene = folderCOSMIC + '/COSMICParsedGene.txt'
    cosmicParsedGene2Mutation = folderCOSMIC + '/COSMICParsedGene2Mutation.txt'
    cosmicParsedMutation = folderCOSMIC + '/COSMICParsedMutation.txt'

    # DrugBank files used.
    DBTargetFasta = folderDB + '/DrugBankApprovedTargets.fasta'
    DBTargetExternalLinks = folderDB + '/ExternalTargetLinks.csv'
    DBXML = folderDB + '/DrugBank.xml'
    DBDrugIDs = folderDB + '/DBDrugs.txt'
    DBTargetIDs = folderDB + '/DBTargets.txt'

    # Ensembl files used.
    ensemblGeneIDs = folderEnsembl + '/EnsemblGeneIDs.txt'
    ensemblExternalIDsOne = folderEnsembl + '/ExternalOne.txt'
    ensemblExternalIDsTwo = folderEnsembl + '/ExternalTwo.txt'  # Two external ID files are needed due to Ensembl Biomart limitations on the number of types of external ID that can be downloaded in one request.
    ensemblTranscripts = folderEnsembl + '/Transcript.txt'
    ensemblParsedTranscripts = folderEnsembl + '/ParsedTranscript.txt'
    ensemblGermSNPResults = folderEnsembl + '/GermSNPs.txt'
    ensemblParsedGermVariants = folderEnsembl + '/ParsedGermVariants.txt'

    # Ensembl Perl API files used.
    ensemblVariationScript = folderEnsemblPerlAPI + '/EnsemblVariant.pl'
    ensemblHomologScript = folderEnsemblPerlAPI + '/EnsemblHomologs.pl'
    ensemblTaxonomyMap = folderEnsembl + '/EnsemblTaxonomyMap.txt'
    ensemblHomologData = folderEnsembl + '/EnsemblHomologs.txt'
    ensemblGenomesHomologScript = folderEnsemblPerlAPI + '/EnsemblGenomesHomologs.pl'
    ensemblGenomesTaxonomyMap = folderEnsembl + '/EnsemblGenomesTaxonomyMap.txt'
    ensemblGenomesHomologData = folderEnsembl + '/EnsemblGenomesHomologs.txt'
    ensemblParsedHomology = folderEnsembl + '/EnsemblParsedHomologs.txt'

    # epestfind files used.
    epestfindExe = folderEpestfind + '/epestfind.exe'

    # Gene Ontology files used.
    completeGODatabase = folderGO + '/go_daily-termdb-data'  # The GO database.
    parsedGOOutput = folderGO + '/GOParsed.txt'  # A file of 6-tuples, one on each line.
                                                 # The first element is the numerical identifier of the GO term.
                                                 # The second element is the name of the GO term.
                                                 # The third element is the category (biological_process, cellular component or molecular_function) that the term belongs to.
                                                 # The fourth element is all the paths from the term to the category it belongs to. The paths are separated from one another using semi-colons, and the elements of each path are separated from one another using '#'.
                                                 # The fifth element is all the level one terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the category in element three.
                                                 # The sixth element is all the level two terms along the paths. These are all the terms that are in a path in element four and are diect descendants of the terms in element five.

    # PathwayCommons files used.
    pathwayElements = folderPathwayCommons + '/PathwayElements.txt'  # The file containing all the pathways in PathwayCommons, and their elements.

    # pepstats files used.
    pepstatsExe = folderPepstats + '/pepstats.exe'

    # SEG files used.
    SEGExe = folderSEG + '/segmasker.exe'

    # TTD files used.
    TTDTargets = folderTTD + '/TTDTargetDataset.txt'  # The TTD drug target database in the raw format.
    TTDUPAccessions = folderTTD + '/UPAccessions.txt'  # A file of UniProt accessions, one on each line.
                                                       # Each UP accession in the file is the target of an approved drug (as recorded by the TTD).
    TTDDrugXref = folderTTD + '/TTDDrugXref.txt'  # The TTD data for the drugs in the database.
    TTDTarget2Drug = folderTTD + '/TTDTarget2Drug.txt'  # A tab separated (tsv) file, with three elements on each line.
                                                        # The first element is the TTD ID for the protein.
                                                        # The second element is a comma separated list of UniProt accessions that the first element is linked to.
                                                        # The third element is a semi-colon separated list of approved drugs that target the protein. For each drug a comma separated list of three elements is recorded.
                                                        #     The first element is the name of the drug (as recorded by the TTD).
                                                        #     The second element is the CAS number of the drug.
                                                        #     The third element is the PubChem CID for the drug.

    # UniGene files used.
    unigeneProfiles = folderUG + '/Hs.profiles'  # The file containing the expression profiles for the human clusters.
    unigeneParsedOutput = folderUG + '/UniGeneParsed.txt'  # A file of tuples, one on each line.
                                                           # The first element of the tuple is the UniGene numerical Id of the cluster.
                                                           # The remaining elements of the tuple are records of the number of ESTs expressed for each expression option.
    unigeneParsedTotals = folderUG + '/UniGeneTotals.txt'  # A file of 2-tuples, one on each line.
                                                           # The first element of the tuple is the body site, developmental stage or health site of the expression.
                                                           # The second element of the tuple is the number of ESTs for the given body site, developmental stage or health state over all cluster.

    # UniProt files used.
    UPHuman = folderUP + '/UniProtHuman.xml'  # The XML file of all UniProt human proteins.
    UPGPCRs = folderUP + '/UniProtGPCRs.txt'  # The file containing all the UniProt accessions of the GPCRs in UniProt.
    UPKinases = folderUP + '/UniProtKinases.txt'  # The file containing the UniProt accessions of the kinases in UniProt.
    UPProteases = folderUP + '/UniProtPeptidases.txt'  # The file containing the UniProt accessions of the proteases in UniProt.
    UPDrugIDs = folderUP + '/DBDrugs.txt'  # A tab separated (tsv) file, with two elements on each line.
                                           # The first element is a representative UniProt accession.
                                           # The second element is a semi-colon separated list of DrugBank drug IDs, one ID for each DrugBank drug that is recorded as being linked to the accession in UniProt.
    UPHumanAccessionMap = folderUP + '/AccMap.txt'  # A tab separated (tsv) file, with two elements on each line.
                                                    # The first element is a UniProt protein accession (representative or not (most often not)).
                                                    # The second elements is the representative UniProt accession that the accession in the first element maps to.
    UPHumanAccessions = folderUP + '/UPAccessions.txt'  # A file with one representative human protein UniProt accession on each line.
                                                        # The number of accessions in this file is the same as the number of human proteins recorded in UniProt.
    UPHumanNames = folderUP + '/UPNames.txt'  # A file with one human protein name on each line.
                                              # The number of names in this file is the same as the number of human proteins recorded in UniProt.
    UPProteinInfo = folderUP + '/UPProteinInfo.txt'  # A file of 56-tuples, one on each line. See the README for a full description of the file.
    UPExternalLinks = folderUP + '/UPExternalLinks.csv'  # A comma separated (csv) file, with six elements on each line.
                                                         # The first element is a representative UniProt accession.
                                                         # The second element is a semi-colon separated list of the Entrez Gene IDs that are recorded as being linked to the accession in UniProt.
                                                         # The third element is a semi-colon separated list of the UniGene cluster IDs that are recorded as being linked to the accession in UniProt.
                                                         # The fourth element is a semi-colon separated list of the Gene Ontology term IDs that are recorded as being linked to the accession in UniProt.
                                                         # The fifth element is a semi-colon separated list of the HGNC IDs that are recorded as being linked to the accession in UniProt.
                                                         # The sixth element is a semi-colon separated list of '-' separated lists of Ensembl records that are recorded as being linked to the accession in UniProt. The format for the '-' separated lists is:
                                                         #     The first element is the Ensembl Gene ID that is recorded as being linked to the UniProt accession (because the transcript in the second element is linked to it).
                                                         #     The second element is the Ensembl Transcript ID that is recorded as being linked to the UniProt accession (because the protein in the third element is linked to it).
                                                         #     The third element is the Ensembl Protein ID that is recorded as being linked to the UniProt accession.
                                                         #     Example : ENSG00000143627-ENST00000342741-ENSP00000339933;ENSG00000143627-ENST00000271946-ENSP00000271946
    UPPPIData = folderUP + '/UPPPIData.txt'  # A file of 5-tuples, one on each line.
                                             # The first element is the UniProt accession for the first protein in the interaction.
                                             # The second element is the UniProt accession for the second protein in the interaction.
                                             # The third element is a number if the second protein is an isofom (the number being the portion after the '-' in a UniProt isoform accession (e.g. 3 in O00257-3)), or 'No Isoform' if the second protein is not an isoform.
                                             # The fourth element is false if the organism the second protein comes from is the same as that of the first (i.e. both are human), or true if the second protein is non-human.
                                             # The fifth element is the number of experiments that give evidence for the interaction between the first and second proteins.

    #===========================================================================
    # Database Schemas and Tables.
    #===========================================================================
    # Names for the MySQL schemas used.
    schemaChEMBL = 'chembl'  # The schema to store the ChEMBL database.
    schemaGO = 'go'  # The schema for the entire GO database.
    schemaPharmGKB = 'pharmgkb'  # The schema for the PharmGKB database.
    schemaProteins = 'proteindatabase'  # The schema used to hold the human protein information tables.

    # Names for the MySQL tables used to store information about the human proteome.
    tableEnsemblGene = schemaProteins + '.ensemblgene'  # The table used to store information about Ensembl genes.
    tableCancerGene = schemaProteins + '.cancergene'  # The table used to store information from Entrez Gene.
    tableGOInfo = schemaProteins + '.goinfo'  # The table used to store information from the GO.
    tableProteinInfo = schemaProteins + '.proteininfo'  # The table used to store information from UniProt.
    tableBLASTResults = schemaProteins + '.blastresults'  # The table used to store the information about the protein similarity.
    tableNonRedundant = schemaProteins + '.nonredundant'  # The table used to store the information about the non-redundant sets of proteins.
    tablePPI = schemaProteins + '.ppi'  # The table used to store the information about the protein protin interactions.
    tableGermVariants = schemaProteins + '.germvariants'  # The table used to store Ensembl information about germ line variants.
    tableHomologs = schemaProteins + '.homologs'  # The table used to store the information about the protein homologs.
    tableUniGene = schemaProteins + '.unigene'  # The table used to store information from UniGene.
    tableUniGeneTotals = schemaProteins + '.unigenetotals'  # The table used to store the totals of the different expression options.
    tableUniProt2Ensembl = schemaProteins + '.uniprot2ensembl'  # The table used to xref UniProt accessions and Ensembl Gene IDs.
    tableUniProt2GO = schemaProteins + '.uniprot2go'  # The table used to xref UniProt accessions and GO term IDs.
    tableUniProt2HGNC = schemaProteins + '.uniprot2hgnc'  # The table used to xref UniProt accessions and HGNC IDs.
    tableUniProt2UniGene = schemaProteins + '.uniprot2unigene'  # The table used to xref UniProt accessions and UniGene IDs.
    tableDrugs = schemaProteins + '.drugs'  # The table used to store the information about the approved drugs.
    tablePathways = schemaProteins + '.pathways'  # The table used to store the information about the pathways.
    tableStability = schemaProteins + '.stability'  # The table used to store the information about the protein stability.
    tableCOSMICGene = schemaProteins + '.cosmicgene'  # The table used to store the COSMIC data.
    tableCOSMICGene2Mutation = schemaProteins + '.cosmicgene2mutation'
    tableCOSMICMutation = schemaProteins + '.cosmicmutation'

    # Names for the MySQL views used.
    viewAllAllTargRN = schemaProteins + '.all_all_targ_r_n'  # A view of all the non-target human proteins.
    viewAllAllTargRP = schemaProteins + '.all_all_targ_r_p'  # A view of all the drug target human proteins.
    viewAllAllTargNRN = schemaProteins + '.all_all_targ_nr_n'  # A view of all the non-redundant non-target human proteins.
    viewAllAllTargNRP = schemaProteins + '.all_all_targ_nr_p'  # A view of all the non-redundant drug target human proteins.

    viewTypeGPCRTargRN = schemaProteins + '.type_gpcr_targ_r_n'  # A view of all the non-target GPCRs.
    viewTypeGPCRTargRP = schemaProteins + '.type_gpcr_targ_r_p'  # A view of all the target GPCRs.
    viewTypeGPCRTargNRN = schemaProteins + '.type_gpcr_targ_nr_n'  # A view of all the non-redundant non-target GPCRs.
    viewTypeGPCRTargNRP = schemaProteins + '.type_gpcr_targ_nr_p'  # A view of all the non-redundant target GPCRs.
    viewTypeIonTargRN = schemaProteins + '.type_ionchannel_targ_r_n'  # A view of all the non-target ion channels.
    viewTypeIonTargRP = schemaProteins + '.type_ionchannel_targ_r_p'  # A view of all the target ion channels.
    viewTypeIonTargNRN = schemaProteins + '.type_ionchannel_targ_nr_n'  # A view of all the non-redundant non-target ion channels.
    viewTypeIonTargNRP = schemaProteins + '.type_ionchannel_targ_nr_p'  # A view of all the non-redundant target ion channels.
    viewTypeKinaseTargRN = schemaProteins + '.type_kinase_targ_r_n'  # A view of all the non-target kinases.
    viewTypeKinaseTargRP = schemaProteins + '.type_kinase_targ_r_p'  # A view of all the target kinases.
    viewTypeKinaseTargNRN = schemaProteins + '.type_kinase_targ_nr_n'  # A view of all the non-redundant non-target kinases.
    viewTypeKinaseTargNRP = schemaProteins + '.type_kinase_targ_nr_p'  # A view of all the non-redundant target kinases.
    viewTypeProteaseTargRN = schemaProteins + '.type_protease_targ_r_n'  # A view of all the non-target proteases.
    viewTypeProteaseTargRP = schemaProteins + '.type_protease_targ_r_p'  # A view of all the target proteases.
    viewTypeProteaseTargNRN = schemaProteins + '.type_protease_targ_nr_n'  # A view of all the non-redundant non-target proteases.
    viewTypeProteaseTargNRP = schemaProteins + '.type_protease_targ_nr_p'  # A view of all the non-redundant target proteases.

    viewIllCancerCTNCNTRN = schemaProteins + '.ill_cancer_ctncnt_r_n'  # A view of all the non-target non-cancer proteins.
    viewIllCancerCTNCNTRP = schemaProteins + '.ill_cancer_ctncnt_r_p'  # A view of all the target cancer proteins.
    viewIllCancerCTNCNTNRN = schemaProteins + '.ill_cancer_ctncnt_nr_n'  # A view of all the non-redundant non-target non-cancer proteins.
    viewIllCancerCTNCNTNRP = schemaProteins + '.ill_cancer_ctncnt_nr_p'  # A view of all the non-redundant target cancer proteins.
    viewIllCancerTargRN = schemaProteins + '.ill_cancer_targ_r_n'  # A view of all the non-target cancer proteins.
    viewIllCancerTargRP = schemaProteins + '.ill_cancer_targ_r_p'  # A view of all the target cancer proteins.
    viewIllCancerTargNRN = schemaProteins + '.ill_cancer_targ_nr_n'  # A view of all the non-redundant non-target cancer proteins.
    viewIllCancerTargNRP = schemaProteins + '.ill_cancer_targ_nr_p'  # A view of all the non-redundant target cancer proteins.
    viewIllCancerTypeRN = schemaProteins + '.ill_cancer_type_r_n'  # A view of all the non-cancer target proteins.
    viewIllCancerTypeRP = schemaProteins + '.ill_cancer_type_r_p'  # A view of all the cancer target proteins.
    viewIllCancerTypeNRN = schemaProteins + '.ill_cancer_type_nr_n'  # A view of all the non-redundant non-cancer target proteins.
    viewIllCancerTypeNRP = schemaProteins + '.ill_cancer_type_nr_p'  # A view of all the non-redundant cancer target proteins.
    viewIllCancerProtRN = schemaProteins + '.ill_cancer_prot_r_n'  # A view of all the non-cancer proteins.
    viewIllCancerProtRP = schemaProteins + '.ill_cancer_prot_r_p'  # A view of all the cancer proteins.
    viewIllCancerProtNRN = schemaProteins + '.ill_cancer_prot_nr_n'  # A view of all the non-redundant non-cancer proteins.
    viewIllCancerProtNRP = schemaProteins + '.ill_cancer_prot_nr_p'  # A view of all the non-redundant cancer proteins.

    #===========================================================================
    # Variables for User Input Parsing.
    #===========================================================================

    # Records whether the database should be initialised.
    doInitialise = False  # doInitialise is True if the user elects to re-initialise the database.
    initialisationScript = DATA + '/Initialisation.sql'

    # Records whether any files/database dumps need loading into the local database.
    # doLoad is True if any files/database dumps should be loaded. toLoad contains the files/database dumps to load.
    doLoad = False  # doLoad is True if files/database dumps should be loaded.
    toLoad = []  # toLoad is a list of the files/database dumps that need loading.
    validLoadArgs = ['GO', 'ChEMBL', 'PharmGKB']  # validLoadArgs is a list of the acronyms of the files/database dumps that can be loaded.

    # Records whether files/tables need parsing, and which ones need parsing.
    doParse = False  # doParse is True if the user selects to parse files/tables.
    toParse = []  # toParse contains the list of files/tables to parse.
    validParseArgs = ['UG', 'CGC', 'GO', 'UP', 'TTD', 'DB', 'BindingDB', 'ChEMBL', 'Ensembl', 'COSMIC', 'CGI']  # validParseArgs is a list of the acronyms of the files/tables that can be parsed.

    # Records whether database tables need updating, and which ones need updating.
    doUpdate = False  # doUpdate is True if the user selects to update database tables.
    toUpdate = []  # toUpdate contains the list of database tables to update.
    validUpdateArgs = ['UG', 'Cancer', 'GO', 'UP', 'Ensembl', 'Xref', 'Target', 'Drug', 'Pathway', 'Stability', 'COSMIC']  # validUpdateArgs is a list of the acronyms of the database tables that can be updates.

    # Records whether the secondary structure prediction need altering.
    doSecondary = False  # doSecondary is True if the user selects to alter secondary structure predictions.
    toSecondary = ''
    outputDirectorySecondary = folderPredictions + '/SecondaryStructure'
    predictionDirectionSecondary = ''
    seqsPerFileSecondary = 100
    maxSeqLenSecondary = 100000

    # Records whether the subcellular localisation predictions need altering.
    doSubcell = False  # doSubcell is True if the user selects to alter subcellular localisation predictions.
    toSubcell = ''  # The table which is to have the predictions altered.
    outputDirectorySubcell = folderPredictions + '/SubcellularLocation'  # The directory from which the predictions should be loaded OR the fasta files written.
    predictionDirectionSubcell = ''  # Either IN (loading predictions into the database), OUTA (all proteins will be predicted) or OUTS (only those proteins that have not been predicted will be output)
    seqsPerFileSubcell = 100  # The number of sequences per fasta file. Only used when the direction is OUTA or OUTS.
    maxSeqLengthSubcell = 100000  # The maximum length of any one sequence in the file. Only used when the direction is OUTA or OUTS.

    # Records whether fasta format files need outputing.
    # The fasta format file of the table toFasta[0] will be saved to the location outputLocationFasta[0].
    doFasta = False  # doFasta is True if the user selects to output fasta files of tables/views.
    toFasta = []  # The list of tables/views to generate fasta files from.
    outputLocationFasta = []  # The list of locations to save the fasta files.

    # Records whether the genetic algorithm data should be generated.
    doGADataGenerate = False  # doGADataGenerate is True if the user selects to generate the GA data.
    viewToGenerateFrom = ''  # The view from which the GA data should be generated from. If you want to get the data about all proteins then you would do all_all_targ
    viewEndings = ['nr_n', 'nr_p']  # These are the endings on the view (i.e. the characters that come after the final '_'. Supplied in the same order as the classification.
    GAClassifications = ['Unlabelled', 'Positive']  # The names for the two classifications.

    #===========================================================================
    # Parse the User Input.
    #===========================================================================
    try:
        flags = 'hil:p:u:s:c:f:g:'
        longFlags = ['--help', '--init', '--load', '--parse', '--update', '--secondary', '--subcell', '--fasta',
                     '--gadata']
        opts, args = getopt.getopt(args, flags, longFlags)
    except getopt.GetoptError:
        print '\nERROR: Unknown command line option found.\n'
        print_help()
        sys.exit()

    for opt, arg in opts:
        if opt in ['-h', '--help']:
            print_help()
            sys.exit()
        elif opt in ['-i', '--init']:
            doInitialise = True
        elif opt in ['-l', '--load']:
            doLoad = True
            toLoad = arg.split('-')
            if not set(toLoad) <= set(validLoadArgs):
                print 'Load ERROR: ', toLoad, validLoadArgs
                sys.exit()
            toLoad = list(set(toLoad))
        elif opt in ['-p', '--parse']:
            doParse = True
            toParse = arg.split('-')
            if not set(toParse) <= set(validParseArgs):
                print 'Parse ERROR: ', toParse, validParseArgs
                sys.exit()
            toParse = list(set(toParse))
        elif opt in ['-u', '--update']:
            doUpdate = True
            toUpdate = arg.split('-')
            if not set(toUpdate) <= set(validUpdateArgs):
                print 'Update ERROR: ', toUpdate, validUpdateArgs
                sys.exit()
            toUpdate = list(set(toUpdate))
        elif opt in ['-s', '--secondary']:
            doSecondary = True
            chunks = arg.split('-')
            toSecondary = chunks[0]
            predictionDirectionSecondary = chunks[1]
            if predictionDirectionSecondary in ['OUTA', 'OUTS']:
                if len(chunks) == 3:
                    seqsPerFileSecondary = int(chunks[2])
                elif len(chunks) == 4:
                    seqsPerFileSecondary = int(chunks[2])
                    maxSeqLenSecondary = int(chunks[3])
                else:
                    print 'Secondary structure prediction ERROR: Not the correct number of arguments'
                    sys.exit()
        elif opt in ['-c',  '--subcell']:
            doSubcell = True
            chunks = arg.split('-')
            toSubcell = chunks[0]
            predictionDirectionSubcell = chunks[1]
            if predictionDirectionSubcell in ['OUTA', 'OUTS']:
                if len(chunks) == 3:
                    seqsPerFileSubcell = int(chunks[2])
                elif len(chunks) == 4:
                    seqsPerFileSubcell = int(chunks[2])
                    maxSeqLengthSubcell = int(chunks[3])
                else:
                    print 'Subcellular location prediction ERROR: Not the correct number of arguments'
        elif opt in ['-f', '--fasta']:
            doFasta = True
            chunks = arg.split('-')
            toFasta = chunks[::2]  # Select all elements of the list with an even index (chunka[0],chunks[2],...).
            outputLocationFasta = chunks[1::2]  # Select all elements of the list with an odd index (chunka[1],chunks[3],...).
        elif opt in ['-g', '--gadata']:
            doGADataGenerate = True
            chunks = arg.split('-')
            viewToGenerateFrom = chunks[0]
            if len(chunks) > 1:
                viewEndings = [chunks[1], chunks[2]]
                GAClassifications = [chunks[3], chunks[4]]

    #===========================================================================
    # Initialising the Database
    #===========================================================================
    if doInitialise:
        confirmation = raw_input('Are you sure you want to re-initialise the database? y for yes and n for no: ')
        if confirmation == 'y' or confirmation == 'Y':
            subprocess.call('mysql.exe -u root -p' + DATABASEPASSWORD + ' mysql < ' + initialisationScript, shell=True, cwd=MYSQLBIN)
            doLoad = True
            toLoad = validLoadArgs
            doParse = True
            toParse = validParseArgs
            doUpdate = True
            toUpdate = validUpdateArgs

    #===========================================================================
    # Loading Database Schemas from Files and Database Dumps.
    #===========================================================================
    if doLoad:
        if 'GO' in toLoad:
            print '\nNow Loading GO'
            loaders.GO.main(DATABASEPASSWORD, schemaGO, completeGODatabase, MYSQLBIN)
        if 'ChEMBL' in toLoad:
            print '\nNow Loading ChEMBL'
            loaders.chembl.main(DATABASEPASSWORD, schemaChEMBL, completeChEMBLDatabase, MYSQLBIN)
        if 'PharmGKB' in toLoad:
            print '\nNow Loading PharmGKB'
            loaders.pharmgkb.main(DATABASEPASSWORD, schemaPharmGKB, pharmGKBDiseases, pharmGKBDrugs, pharmGKBGenes,
                                  pharmGKBRelationships, pharmGKBInitialisation, MYSQLBIN)

    #===========================================================================
    # Parsing Files and Database Tables.
    #===========================================================================
    if doParse:
        if 'UG' in toParse:
            print '\nNow Parsing UniGene'
            parsers.parseUG.main(unigeneProfiles, unigeneParsedOutput, unigeneParsedTotals)
        if 'CGC' in toParse:
            print '\nNow Parsing Cancer Gene Census'
            parsers.parseCGC.main(CGCData, CGCParsed)
        if 'GO' in toParse:
            print '\nNow Parsing GO'
            # Process the GO database in order to extract the information about the paths of each term to its root term.
            parsers.parseGO.main(schemaGO, parsedGOOutput)
        if 'UP' in toParse:
            print '\nNow Parsing Uniprot'
            # Parse the UniProt human proteome file.
            parsers.parseUP.main(UPHuman, UPGPCRs, UPKinases, UPProteases, UPHumanAccessions, UPHumanAccessionMap,
                                 UPHumanNames, UPDrugIDs, UPProteinInfo, UPExternalLinks, UPPPIData, folderUP)
            print '\tNow Updating Cross-references'
            # Extract xrefs from Ensembl.
            martName = 'ensembl'
            datasetName = 'hsapiens_gene_ensembl'
            filterName = 'uniprot_swissprot_accession'
            filterValue = ','.join(utilities.file2list.main(UPHumanAccessions))
            attributes = ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'entrezgene',
                          'unigene', 'uniprot_swissprot_accession']
            utilities.biomartquery.single_filter(martName, datasetName, filterName, filterValue, attributes,
                                                 biomartScript, biomartQuery, ensemblExternalIDsOne)
            attributes = ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',
                          'uniprot_swissprot_accession', 'hgnc_id']
            utilities.biomartquery.single_filter(martName, datasetName, filterName, filterValue, attributes,
                                                 biomartScript, biomartQuery, ensemblExternalIDsTwo)
            parsers.parseEnsembl.update_xref_and_ensembl_IDs(ensemblExternalIDsOne, ensemblExternalIDsTwo, UPExternalLinks, ensemblGeneIDs)
        if 'TTD' in toParse:
            print '\nNow Parsing TTD'
            # Parse the TTD database, to find approved targets.
            parsers.parseTTD.main(TTDTargets, TTDUPAccessions, TTDDrugXref, TTDTarget2Drug)
        if 'DB' in toParse:
            print '\nNow Parsing DrugBank'
            # Parse the DrugBank files, to find the approved targets and drugs.
            parsers.parseDB.main(DBTargetFasta, DBXML, DBTargetExternalLinks, DBTargetIDs, DBDrugIDs)
        if 'BindingDB' in toParse:
            print '\nNow Parsing BindingDB'
            # Parse the BindingDB SDF file to get binding affinity information about potential targets.
            parsers.parseBindingDB.main(bindingSDF, binding2PubChem, bindingParsed)
        if 'ChEMBL' in toParse:
            print '\nNow Parsing ChEMBL'
            # Parse the ChEMBL database, to find approved targets.
            parsers.parseChEMBL.main(ChEMBLUPAccessions, ChEMBLCID, DATABASEPASSWORD, schemaChEMBL)
        if 'Ensembl' in toParse:
            print '\nNow Parsing Ensembl'
            martName = 'ensembl'
            ensemblGeneIDList = sorted(utilities.file2list.main(ensemblGeneIDs))
            ensemblGenes = ','.join(ensemblGeneIDList)

            # Extract the transcript information from Ensembl.
            datasetName = 'hsapiens_gene_ensembl'
            filterName = 'ensembl_gene_id'
            filterValue = ensemblGenes
            attributes = ['ensembl_gene_id', 'ensembl_transcript_id', 'transcript_count', 'transcript_biotype']
            utilities.biomartquery.single_filter(martName, datasetName, filterName, filterValue, attributes,
                                                 biomartScript, biomartQuery, ensemblTranscripts)

            # Extract the germ variant information from Ensembl.
            subprocess.call(['perl', ensemblVariationScript, ensemblGeneIDs, ensemblGermSNPResults])

            # Parse the information from Ensembl.
            parsers.parseEnsembl.main(ensemblTranscripts, ensemblParsedTranscripts, ensemblGermSNPResults, ensemblParsedGermVariants)

            # Extract and parse the homolog information from Ensembl.
           parsers.parseHomologs.main(ensemblHomologScript, ensemblTaxonomyMap, ensemblHomologData,
                                       ensemblGenomesHomologScript, ensemblGenomesTaxonomyMap,
                                       ensemblGenomesHomologData, ensemblParsedHomology, ensemblGeneIDs)
        if 'COSMIC' in toParse:
            print '\nNow Parsing COSMIC'
            parsers.parseCOSMIC.main(cosmicData, cosmicParsedGene, cosmicParsedMutation, cosmicParsedGene2Mutation)
        if 'CGI' in toParse:
            print '\nNow Parsing Cancer Gene Index'
            parsers.parseCGI.main(CGIData, CGIHGNCIDs, CGIUPAccessions)

    #===========================================================================
    # Update the Information in the Database.
    #===========================================================================
    if doUpdate:
        if 'UG' in toUpdate:
            print '\nNow Updating UniGene'
            updaters.updateUG.main(unigeneParsedOutput, unigeneParsedTotals, schemaProteins, tableUniGene,
                                   tableUniGeneTotals, DATABASEPASSWORD)
        if 'COSMIC' in toUpdate:
            print '\nNow Updating COSMIC Information'
            updaters.updateCOSMIC.main(cosmicParsedGene, cosmicParsedGene2Mutation, cosmicParsedMutation,
                                       DATABASEPASSWORD, schemaProteins, tableCOSMICGene, tableCOSMICGene2Mutation, tableCOSMICMutation)
        if 'Cancer' in toUpdate:
            print '\nNow Updating Cancer Information'
            updaters.updateCancer.main(CGCParsed, cancerTargets, CGIHGNCIDs, CGIUPAccessions, UPExternalLinks, UPHumanAccessionMap,
                                       schemaProteins, tableCancerGene, tableCOSMICGene, tableCOSMICGene2Mutation, tableCOSMICMutation,
                                       DATABASEPASSWORD)
        if 'GO' in toUpdate:
            print '\nNow Updating GO'
            updaters.updateGO.main(parsedGOOutput, schemaProteins, tableGOInfo, DATABASEPASSWORD)
        if 'UP' in toUpdate:
            print '\nNow Updating Uniprot'
            updaters.updateUP.main(UPProteinInfo, schemaProteins, tableProteinInfo, folderSEG, SEGExe, folderEpestfind,
                                   epestfindExe, folderPepstats, pepstatsExe, folderBLAST, psiblastExe,
                                   makeBLASTDatabaseExe, tableBLASTResults, DATABASEPASSWORD)
            print '\nNow Updating Protein Protein Interactions'
            updaters.updatePPI.main(UPPPIData, schemaProteins, tablePPI, DATABASEPASSWORD)
        if 'Ensembl' in toUpdate:
            print '\nNow Updating Ensembl'
            updaters.updateEnsembl.main(ensemblParsedTranscripts, ensemblParsedGermVariants, ensemblParsedHomology, schemaProteins,
                                        tableEnsemblGene, tableGermVariants, tableHomologs, DATABASEPASSWORD)
        if 'Xref' in toUpdate:
            print '\nNow Updating Cross-reference Tables'
            # Create a dictionary of the tables from which the IDs of the different databases need to be extracted.
            # In order to update the pivot tables it is necessary to ensure that the UniProt ID is being linked to
            # a valid ID from another database. The indices for the dictionary are the columns from which to extract
            # the information. The entry for each index is the table from which to extract the index column.
            tableDict = {'EnsemblGeneID' : [tableEnsemblGene], 'GOTermID' : [tableGOInfo],
                         'UPAccession' : [tableProteinInfo], 'UniGeneID' : [tableUniGene]}
            updaters.updatexref.main(UPExternalLinks, tableUniProt2Ensembl, tableUniProt2GO,
                                     tableUniProt2UniGene, tableUniProt2HGNC, tableDict, schemaProteins, DATABASEPASSWORD)
        if 'Target' in toUpdate:
            print '\nNow Updating TTD, ChEMBL, DrugBank and Redundancy Information.'
            # Create a dictionary of the views that store the UP accession and sequence information for the proteins
            # of different target/non-target classifications. The dictionary is indexed by the names of the columns in
            # the nonredundant table. The value stored at an index is the view which holds the proteins that need to
            # be made non-redundant for the particular index column.
            viewsDict = {'AllTargetNegative' : viewAllAllTargRN, 'AllTargetPositive' : viewAllAllTargRP,
                         'GPCRTargetNegative' : viewTypeGPCRTargRN, 'GPCRTargetPositive' : viewTypeGPCRTargRP,
                         'IonChannelTargetNegative' : viewTypeIonTargRN, 'IonChannelTargetPositive' : viewTypeIonTargRP,
                         'KinaseTargetNegative' : viewTypeKinaseTargRN, 'KinaseTargetPositive' : viewTypeKinaseTargRP,
                         'ProteaseTargetNegative' : viewTypeProteaseTargRN, 'ProteaseTargetPositive' : viewTypeProteaseTargRP,
                         'CancerCTNCNTNegative' : viewIllCancerCTNCNTRN, 'CancerCTNCNTPositive' : viewIllCancerCTNCNTRP,
                         'CancerTargetNegative' : viewIllCancerTargRN, 'CancerTargetPositive' : viewIllCancerTargRP,
                         'CancerTypeNegative' : viewIllCancerTypeRN, 'CancerTypePositive' : viewIllCancerTypeRP,
                         'CancerProteinNegative' : viewIllCancerProtRN, 'CancerProteinPositive' : viewIllCancerProtRP
                         }
            updaters.updatetargetandredundancy.main(DBDrugIDs, DBTargetIDs, TTDUPAccessions,
                                                    ChEMBLUPAccessions, UPHumanAccessionMap, UPDrugIDs, pharmGKBDrugTargets, folderCulling,
                                                    schemaProteins, tableProteinInfo, tableNonRedundant,
                                                    tableBLASTResults, DATABASEPASSWORD, viewsDict)
        if 'Drug' in toUpdate:
            print '\nNow Updating Drug Information'
            updaters.updatedrug.main(UPDrugIDs, DBDrugIDs, DBTargetIDs, TTDTarget2Drug, ChEMBLUPAccessions, ChEMBLCID,
                                     bindingParsed, UPHumanAccessionMap, DATABASEPASSWORD, schemaProteins, tableDrugs)
        if 'Pathway' in toUpdate:
            print '\nNow Updating Pathway Information'
            updaters.updatepathway.main(pathwayElements, UPHumanNames, UPHumanAccessionMap, DATABASEPASSWORD,
                                        schemaProteins, tablePathways)
        if 'Stability' in toUpdate:
            print '\nNow Updating Protein Stability Information'
            updaters.updatestability.main(DATABASEPASSWORD, schemaProteins, tableProteinInfo, tableStability)

    #===========================================================================
    # Alter Secondary Structure Prediction Information.
    #===========================================================================
    if doSecondary:
        utilities.predictions.main(toSecondary, outputDirectorySecondary, predictionDirectionSecondary,
                                   seqsPerFileSecondary, maxSeqLenSecondary, schemaProteins, DATABASEPASSWORD,
                                   'PredictedSecondaryStructure')

    #===========================================================================
    # Alter Subcellular Localisation Prediction Information.
    #===========================================================================
    if doSubcell:
        utilities.predictions.main(toSubcell, outputDirectorySubcell, predictionDirectionSubcell, seqsPerFileSubcell,
                                   maxSeqLengthSubcell, schemaProteins, DATABASEPASSWORD,
                                   'PredictedSubcellularLocation')

    #===========================================================================
    # Output Fasta Files of the Specified Tables and/or Views.
    #===========================================================================
    if doFasta:
        print '\nCreating Fasta Files.'
        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        for i in range(len(toFasta)):
            print '\tNow making fasta file from table: ', toFasta[i]
            try:
                cursor = mysql.tableSELECT(cursor, 'UPAccession, Sequence', toFasta[i])
                results = cursor.fetchall()
                utilities.list2file.main(['>' + '\n'.join([j[0], j[1]]) for j in results], outputLocationFasta[i])
            except:
                print '\t\tERROR: Creation of fasta file for table', toFasta[i], ' failed. Please check that the table',
                print '\t\tcontains columns titled UPAccession and Sequence, and that the file location is valid.'
        mysql.closeConnection(conn, cursor)

    #===========================================================================
    # Generate the Data File Needed for the Running of the Genetic Algorithm
    #===========================================================================
    if doGADataGenerate:
        # Assumes that the columns are identical for both tables/views.
        print '\nGenerating GA Data File.'
        outputDirectory = folderGAData + '/' + viewToGenerateFrom.upper()
        if os.path.isdir(outputDirectory):
            shutil.rmtree(outputDirectory)
        os.mkdir(outputDirectory)
        geneticAlgorithmProcessedData = outputDirectory + '/' + viewToGenerateFrom.upper() + '.txt'
        categoricalMappingDataLocation = outputDirectory + '/CategoricalMapping.txt'
        columnDataLocation = outputDirectory + '/Columns.txt'
        ECDataLocation = outputDirectory + '/ECNumbers.txt'
        subcellLocation = outputDirectory + '/SubcellularLocation.txt'
        healthStateLocation = outputDirectory + '/HealthState.txt'
        bodySiteLocation = outputDirectory + '/BodySite.txt'
        developmentalStageLocation = outputDirectory + '/DevelopmentalStage.txt'

        viewDict = {viewAllAllTargNRN : 'AllTargetNegative', viewAllAllTargNRP : 'AllTargetPositive',
                    viewTypeGPCRTargNRN : 'GPCRTargetNegative', viewTypeGPCRTargNRP : 'GPCRTargetPositive',
                    viewTypeIonTargNRN : 'IonChannelTargetNegative', viewTypeIonTargNRP : 'IonChannelTargetPositive',
                    viewTypeKinaseTargNRN : 'KinaseTargetNegative', viewTypeKinaseTargNRP : 'KinaseTargetPositive',
                    viewTypeProteaseTargNRN : 'ProteaseTargetNegative', viewTypeProteaseTargNRP : 'ProteaseTargetPositive',
                    viewIllCancerCTNCNTNRN : 'CancerCTNCNTNegative', viewIllCancerCTNCNTNRP : 'CancerCTNCNTPositive',
                    viewIllCancerTargNRN : 'CancerTargetNegative', viewIllCancerTargNRP : 'CancerTargetPositive',
                    viewIllCancerTypeNRN : 'CancerTypeNegative', viewIllCancerTypeNRP : 'CancerTypePositive',
                    viewIllCancerProtNRN : 'CancerProteinNegative', viewIllCancerProtNRP : 'CancerProteinPositive'
                    }
        positiveColumn = viewDict[schemaProteins + '.' + viewToGenerateFrom + '_' + viewEndings[1]]
        negativeColumn = viewDict[schemaProteins + '.' + viewToGenerateFrom + '_' + viewEndings[0]]

        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        cursor.execute('SHOW COLUMNS FROM ' + viewToGenerateFrom + '_' + viewEndings[1])
        columns = cursor.fetchall()
        columns = [i[0] for i in columns]

        nonredColumns = cursor.execute('SHOW COLUMNS FROM ' + tableNonRedundant)
        nonredColumns = cursor.fetchall()
        nonredColumns = [i[0] for i in nonredColumns]
        nonredundantProteins = cursor.execute('SELECT * FROM ' + tableNonRedundant)
        nonredundantProteins = cursor.fetchall()
        positiveProteinAccs = [i[0] for i in nonredundantProteins if i[nonredColumns.index(positiveColumn)] == 'Y']
        negativeProteinAccs = [i[0] for i in nonredundantProteins if i[nonredColumns.index(negativeColumn)] == 'Y']

        protColumns = cursor.execute('SHOW COLUMNS FROM ' + tableProteinInfo)
        protColumns = cursor.fetchall()
        protColumns = [i[0] for i in protColumns]
        proteins = cursor.execute('SELECT * FROM ' + tableProteinInfo)
        proteins = cursor.fetchall()
        positiveProteins = dict([(i[0], i) for i in proteins if i[0] in positiveProteinAccs])
        negativeProteins = dict([(i[0], i) for i in proteins if i[0] in negativeProteinAccs])

        stabilityColumns = cursor.execute('SHOW COLUMNS FROM ' + tableStability)
        stabilityColumns = cursor.fetchall()
        stabilityColumns = [i[0] for i in stabilityColumns]
        stability = cursor.execute('SELECT * FROM ' + tableStability)
        stability = cursor.fetchall()
        positiveStability = dict([(i[0], i) for i in stability if i[0] in positiveProteinAccs])
        negativeStability = dict([(i[0], i) for i in stability if i[0] in negativeProteinAccs])

        ppi = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_ppi')
        ppi = cursor.fetchall()
        positivePPI = dict([(i[0], i) for i in ppi if i[0] in positiveProteinAccs])
        negativePPI = dict([(i[0], i) for i in ppi if i[0] in negativeProteinAccs])

        expressionColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_expression')
        expressionColumns = cursor.fetchall()
        expressionColumns = [i[0] for i in expressionColumns]
        expression = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_expression')
        expression = cursor.fetchall()
        positiveExpression = dict([(i[0], i) for i in expression if i[0] in positiveProteinAccs])
        negativeExpression = dict([(i[0], i) for i in expression if i[0] in negativeProteinAccs])

        paralogColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_paralogs')
        paralogColumns = cursor.fetchall()
        paralogColumns = [i[0] for i in paralogColumns]
        paralog = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_paralogs')
        paralog = cursor.fetchall()
        positiveParalog = dict([(i[0], i) for i in paralog if i[0] in positiveProteinAccs])
        negativeParalog = dict([(i[0], i) for i in paralog if i[0] in negativeProteinAccs])

        transcriptColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_transcripts')
        transcriptColumns = cursor.fetchall()
        transcriptColumns = [i[0] for i in transcriptColumns]
        transcript = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_transcripts')
        transcript = cursor.fetchall()
        positiveTranscript = dict([(i[0], i) for i in transcript if i[0] in positiveProteinAccs])
        negativeTranscript = dict([(i[0], i) for i in transcript if i[0] in negativeProteinAccs])

        variantColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_germvariants')
        variantColumns = cursor.fetchall()
        variantColumns = [i[0] for i in paralogColumns]
        variant = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_germvariants')
        variant = cursor.fetchall()
        positiveVariant = {}
        negativeVariant = {}
        for i in variant:
            if i[0] in positiveProteinAccs:
                if positiveVariant.has_key(i[0]):
                    positiveVariant[i[0]]['3untrans'] += i[2]
                    positiveVariant[i[0]]['5untrans'] += i[3]
                    positiveVariant[i[0]]['nonsynon'] += i[4]
                    positiveVariant[i[0]]['synon'] += i[5]
                else:
                    positiveVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}
            elif i[0] in negativeProteinAccs:
                if negativeVariant.has_key(i[0]):
                    negativeVariant[i[0]]['3untrans'] += i[2]
                    negativeVariant[i[0]]['5untrans'] += i[3]
                    negativeVariant[i[0]]['nonsynon'] += i[4]
                    negativeVariant[i[0]]['synon'] += i[5]
                else:
                    negativeVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}

        mysql.closeConnection(conn, cursor)

        resultsTarget = []
        resultsNonTarget = []
        for i in [[positiveProteinAccs, 'positive'], [negativeProteinAccs, 'negative']]:
            for j in i[0]:
                currentRecord = []

                # Protein properties.
                if i[1] == 'positive':
                    for k in ['UPAccession', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                              'NegativelyCharged', 'PositivelyCharged', 'Basic', 'Charged', 'Polar', 'NonPolar', 'Aromatic', 'Aliphatic', 'Small', 'Tiny',
                              'PESTMotif', 'LowComplexity', 'Hydrophobicity', 'Isoelectric', 'ECNumber', 'OGlycosylation', 'NGlycosylation', 'Phosphoserine',
                              'Phosphothreonine', 'Phosphotyrosine', 'SubcellularLocation', 'TopologicalDomain', 'PredictedSubcellularLocation',
                              'SignalPeptide', 'TransmembraneHelices', 'AlphaHelices', 'BetaStrands', 'PredictedAlphaHelices', 'PredictedBetaSheets', 'Sequence']:
                        currentRecord.append(positiveProteins[j][protColumns.index(k)])
                else:
                    for k in ['UPAccession', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                              'NegativelyCharged', 'PositivelyCharged', 'Basic', 'Charged', 'Polar', 'NonPolar', 'Aromatic', 'Aliphatic', 'Small', 'Tiny',
                              'PESTMotif', 'LowComplexity', 'Hydrophobicity', 'Isoelectric', 'ECNumber', 'OGlycosylation', 'NGlycosylation', 'Phosphoserine',
                              'Phosphothreonine', 'Phosphotyrosine', 'SubcellularLocation', 'TopologicalDomain', 'PredictedSubcellularLocation',
                              'SignalPeptide', 'TransmembraneHelices', 'AlphaHelices', 'BetaStrands', 'PredictedAlphaHelices', 'PredictedBetaSheets', 'Sequence']:
                        currentRecord.append(negativeProteins[j][protColumns.index(k)])

                # Expression properties.
                if i[1] == 'positive':
                    for k in ['DS_Embryoid_Body', 'DS_Blastocyst', 'DS_Fetus', 'DS_Neonate', 'DS_Infant', 'DS_Juvenile', 'DS_Adult', 'HS_Adrenal_Tumor',
                              'HS_Bladder_Carcinoma', 'HS_Breast_Mammary_Gland_Tumor', 'HS_Cervical_Tumor', 'HS_Chondrosarcoma', 'HS_Colorectal_Tumor',
                              'HS_Esophageal_Tumor', 'HS_Gastrointestinal_Tumor', 'HS_Germ_Cell_Tumor', 'HS_Glioma', 'HS_Head_And_Neck_Tumor',
                              'HS_Kidney_Tumor', 'HS_Leukemia_Tumor', 'HS_Liver_Tumor', 'HS_Lung_Tumor', 'HS_Lymphoma', 'HS_Non_neoplasia', 'HS_Normal',
                              'HS_Ovarian_Tumor', 'HS_Pancreatic_Tumor', 'HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS', 'HS_Prostate_Cancer',
                              'HS_Retinoblastoma', 'HS_Skin_Tumor', 'HS_Soft_Tissue_Muscle_Tissue_Tumor', 'HS_Uterine_Tumor', 'BS_Adipose_Tissue',
                              'BS_Adrenal_Gland', 'BS_Ascites', 'BS_Bladder', 'BS_Blood', 'BS_Bone', 'BS_Bone_Marrow', 'BS_Brain', 'BS_Cervix',
                              'BS_Connective_Tissue', 'BS_Ear', 'BS_Embryonic_Tissue', 'BS_Esophagus', 'BS_Eye', 'BS_Heart', 'BS_Intestine', 'BS_Kidney',
                              'BS_Larynx', 'BS_Liver', 'BS_Lung', 'BS_Lymph', 'BS_Lymph_Node', 'BS_Mammary_Gland', 'BS_Mouth', 'BS_Muscle', 'BS_Nerve',
                              'BS_Ovary', 'BS_Pancreas', 'BS_Parathyroid', 'BS_Pharynx', 'BS_Pituitary_Gland', 'BS_Placenta', 'BS_Prostate',
                              'BS_Salivary_Gland', 'BS_Skin', 'BS_Spleen', 'BS_Stomach', 'BS_Testis', 'BS_Thymus', 'BS_Thyroid', 'BS_Tonsil', 'BS_Trachea',
                              'BS_Umbilical_Cord', 'BS_Uterus', 'BS_Vascular']:
                        currentRecord.append(int(positiveExpression[j][expressionColumns.index(k)]))
                else:
                    for k in ['DS_Embryoid_Body', 'DS_Blastocyst', 'DS_Fetus', 'DS_Neonate', 'DS_Infant', 'DS_Juvenile', 'DS_Adult', 'HS_Adrenal_Tumor',
                              'HS_Bladder_Carcinoma', 'HS_Breast_Mammary_Gland_Tumor', 'HS_Cervical_Tumor', 'HS_Chondrosarcoma', 'HS_Colorectal_Tumor',
                              'HS_Esophageal_Tumor', 'HS_Gastrointestinal_Tumor', 'HS_Germ_Cell_Tumor', 'HS_Glioma', 'HS_Head_And_Neck_Tumor',
                              'HS_Kidney_Tumor', 'HS_Leukemia_Tumor', 'HS_Liver_Tumor', 'HS_Lung_Tumor', 'HS_Lymphoma', 'HS_Non_neoplasia', 'HS_Normal',
                              'HS_Ovarian_Tumor', 'HS_Pancreatic_Tumor', 'HS_Primitive_Neuroectodermal_Tumor_Of_The_CNS', 'HS_Prostate_Cancer',
                              'HS_Retinoblastoma', 'HS_Skin_Tumor', 'HS_Soft_Tissue_Muscle_Tissue_Tumor', 'HS_Uterine_Tumor', 'BS_Adipose_Tissue',
                              'BS_Adrenal_Gland', 'BS_Ascites', 'BS_Bladder', 'BS_Blood', 'BS_Bone', 'BS_Bone_Marrow', 'BS_Brain', 'BS_Cervix',
                              'BS_Connective_Tissue', 'BS_Ear', 'BS_Embryonic_Tissue', 'BS_Esophagus', 'BS_Eye', 'BS_Heart', 'BS_Intestine', 'BS_Kidney',
                              'BS_Larynx', 'BS_Liver', 'BS_Lung', 'BS_Lymph', 'BS_Lymph_Node', 'BS_Mammary_Gland', 'BS_Mouth', 'BS_Muscle', 'BS_Nerve',
                              'BS_Ovary', 'BS_Pancreas', 'BS_Parathyroid', 'BS_Pharynx', 'BS_Pituitary_Gland', 'BS_Placenta', 'BS_Prostate',
                              'BS_Salivary_Gland', 'BS_Skin', 'BS_Spleen', 'BS_Stomach', 'BS_Testis', 'BS_Thymus', 'BS_Thyroid', 'BS_Tonsil', 'BS_Trachea',
                              'BS_Umbilical_Cord', 'BS_Uterus', 'BS_Vascular']:
                        currentRecord.append(int(negativeExpression[j][expressionColumns.index(k)]))

                # SNPs.
                if i[1] == 'positive':
                    if positiveVariant.has_key(j):
                        currentRecord.append(positiveVariant[j]['3untrans'])
                        currentRecord.append(positiveVariant[j]['5untrans'])
                        currentRecord.append(positiveVariant[j]['nonsynon'])
                        currentRecord.append(positiveVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                else:
                    if negativeVariant.has_key(j):
                        currentRecord.append(negativeVariant[j]['3untrans'])
                        currentRecord.append(negativeVariant[j]['5untrans'])
                        currentRecord.append(negativeVariant[j]['nonsynon'])
                        currentRecord.append(negativeVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)

                # Paralogs.
                if i[1] == 'positive':
                    if positiveParalog.has_key(j):
                        currentRecord.append(int(positiveParalog[j][-1]))
                    else:
                        currentRecord.append(0)
                else:
                    if negativeParalog.has_key(j):
                        currentRecord.append(int(negativeParalog[j][-1]))
                    else:
                        currentRecord.append(0)

                # PPI.
                if i[1] == 'positive':
                    if positivePPI.has_key(j):
                        currentRecord.append(int(positivePPI[j][-1]))
                    else:
                        currentRecord.append(0)
                else:
                    if negativePPI.has_key(j):
                        currentRecord.append(int(negativePPI[j][-1]))
                    else:
                        currentRecord.append(0)

                # Transcripts.
                if i[1] == 'positive':
                    if positiveTranscript.has_key(j):
                        currentRecord.append(int(positiveTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)
                else:
                    if negativeTranscript.has_key(j):
                        currentRecord.append(int(negativeTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)

                # Stability properties.
                if i[1] == 'positive':
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(positiveStability[j][stabilityColumns.index(k)])
                else:
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(negativeStability[j][stabilityColumns.index(k)])

                if i[1] == 'positive':
                    resultsTarget.append(tuple(currentRecord))
                else:
                    resultsNonTarget.append(tuple(currentRecord))

#        utilities.gadatageneration.fortran(resultsTarget, resultsNonTarget, columns, geneticAlgorithmProcessedData, columnDataLocation, ECDataLocation,
#                                           subcellLocation, healthStateLocation, bodySiteLocation, developmentalStageLocation)
#        utilities.gadatageneration.fortran_split(geneticAlgorithmProcessedData, classes=[1, 2], splits=5, outputLocation=outputDirectory)
        utilities.gadatageneration.pulearning(resultsTarget, resultsNonTarget, columns, geneticAlgorithmProcessedData, categoricalMappingDataLocation,
                                           columnDataLocation, ECDataLocation, subcellLocation, healthStateLocation, bodySiteLocation,
                                           developmentalStageLocation)

        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        # Get the GO term data for the target proteins.
        targetProteinAccs = [i[0] for i in resultsTarget]
        cursor = mysql.tableSELECT(cursor, '*', tableUniProt2GO, 'UPAccession IN ("' + '","'.join(targetProteinAccs) + '")')
        results = cursor.fetchall()
        targetGOIDs = {}
        for i in results:
            UPAccession = i[0]
            GOTermID = i[1]
            if targetGOIDs.has_key(UPAccession):
                targetGOIDs[UPAccession].add(GOTermID)
            else:
                targetGOIDs[UPAccession] = set([GOTermID])
        targetGOTerms = {}
        for i in targetGOIDs.keys():
            UPAccession = i
            targetGOTerms[UPAccession] = {'biological_process' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                          'cellular_component' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                          'molecular_function' : {'LevelOne' : set([]), 'LevelTwo' : set([])}
                                          }
            GOIDs = [str(i) for i in targetGOIDs[i]]
            cursor = mysql.tableSELECT(cursor, 'GOType, LevelOne, LevelTwo', tableGOInfo, 'GOTermID IN ("' + '","'.join(GOIDs) + '")')
            results = cursor.fetchall()
            for j in results:
                GOType = j[0]
                levelOne = j[1].split(';')
                levelTwo = j[2].split(';')

                if levelOne != ['NA']:
                    targetGOTerms[UPAccession][GOType]['LevelOne'].update(levelOne)
                if levelTwo != ['NA']:
                    targetGOTerms[UPAccession][GOType]['LevelTwo'].update(levelTwo)

        # Get the GO term data for the non-target proteins.
        nonTargetProteinAccs = [i[0] for i in resultsNonTarget]
        cursor = mysql.tableSELECT(cursor, '*', tableUniProt2GO, 'UPAccession IN ("' + '","'.join(nonTargetProteinAccs) + '")')
        results = cursor.fetchall()
        nonTargetGOIDs = {}
        for i in results:
            UPAccession = i[0]
            GOTermID = i[1]
            if nonTargetGOIDs.has_key(UPAccession):
                nonTargetGOIDs[UPAccession].add(GOTermID)
            else:
                nonTargetGOIDs[UPAccession] = set([GOTermID])
        nonTargetGOTerms = {}
        for i in nonTargetGOIDs.keys():
            UPAccession = i
            nonTargetGOTerms[UPAccession] = {'biological_process' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                             'cellular_component' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                             'molecular_function' : {'LevelOne' : set([]), 'LevelTwo' : set([])}
                                             }
            GOIDs = [str(i) for i in nonTargetGOIDs[i]]
            cursor = mysql.tableSELECT(cursor, 'GOType, LevelOne, LevelTwo', tableGOInfo, 'GOTermID IN ("' + '","'.join(GOIDs) + '")')
            results = cursor.fetchall()
            for j in results:
                GOType = j[0]
                levelOne = j[1].split(';')
                levelTwo = j[2].split(';')

                if levelOne != ['NA']:
                    nonTargetGOTerms[UPAccession][GOType]['LevelOne'].update(levelOne)
                if levelTwo != ['NA']:
                    nonTargetGOTerms[UPAccession][GOType]['LevelTwo'].update(levelTwo)

        mysql.closeConnection(conn, cursor)

        utilities.generateGOsummarydata.main(targetGOTerms, nonTargetGOTerms, outputDirectory)

        subprocess.call(['Rscript.exe', folderGAData + '/StatisticalTests.R', outputDirectory, GAClassifications[0], GAClassifications[1]])

if __name__ == '__main__':
    main(sys.argv[1:])