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
    folderBiomart = folderEnsembl + '/Biomart'
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
    bindingSDF = folderBindingDB + '/BindingDB2D.sdf'  # A file containing information about interactions between compounds and targets.
    binding2PubChem = folderBindingDB + '/BDB_cid.txt'  # A file containing mappings from BindingDb compound IDs to PubChem CIDs.
    bindingParsed = folderBindingDB + '/BindingDBParsed.txt'  # A tab separated file (tsv), with four elements on each line.
                                                              # The first element is the UniProt accessions targeted by a compound.
                                                              # The second element is the CID of the compound that targets the UniProt accessions in the first element.
                                                              # The third element is the Ki for the interaction between the compound in element two and the proteins in element one.
                                                              # The fourth element is the Kd for the interaction between the compound in element two and the proteins in element one.

    # Biomart files used.
    biomartScript = folderBiomart + '/biomartquery.pl'  # The location of the script that submits queries to Biomart.
    biomartQuery = folderBiomart + '/BiomartQuery.xml'  # The location where all Biomart XML queries are stored before being submitted.

    # BLAST files used.
    psiblastExe = folderBLAST + '/psiblast.exe'  # Location of the PSI-BLAST executable.
    makeBLASTDatabaseExe = folderBLAST + '/makeblastdb.exe'  # Location of the executable to make a BLAST database.

    # File of all known FDA approved cancer drug targets.
    cancerTargets = folderCancerTargets + '/CancerDrugs.txt'  # Location of the file containing all the known FDA approved cancer targets.

    # Cancer Gene Census files used.
    CGCData = folderCGC + '/cancer_gene_census.tsv'  # A file containing the Cancer Gene Census data.
    CGCParsed = folderCGC + '/CGCParsedData.txt'  # A tab separated (tsv) file, with three elements on each line.
                                                  # The first element is the NCBI gene ID of the gene.
                                                  # The second element is Y if the gene contains a cancer causing germline mutation, else it is N.
                                                  # The third element is Y if the gene contains a cancer causing somatic mutation, else it is N.

    # Cancer Gene Index files used.
    CGIData = folderCGI + '/NCI_CancerIndex_allphases_disease.xml'  # A file containing the Cancer Gene Index data.
    CGIUPAccessions = folderCGI + '/CGIUPAccessions.txt'  # A file with one UniProt accession on each line.
    CGIHGNCIDs = folderCGI + '/CGIHGNCIDs.txt'  # A file with one HGNC gene ID on each line.

    # ChEMBL files used.
    completeChEMBLDatabase = folderChEMBL + '/ChEMBL.sql'  # The ChEMBL database.
    ChEMBLUPAccessions = folderChEMBL + '/UPAccessions.txt'  # A tab separated (tsv) file, with seven elements on each line. There is one line for each protein that an approved target targets.
                                                             # The first element is the UniProt accession for the target protein.
                                                             # The second element is ChEMBL ID of the compound.
                                                             # The third element is the name of the protein in the first element.
                                                             # The fourth element is the activity relation (=, <, <=, > or >=).
                                                             # The fifth element is the value of the activity observed between the compound and the protein.
                                                             # The sixth element is the units of the value in the fifth element.
                                                             # The seventh element is the type of activity measured (Ki, Kd, EC50, etc.)
    ChEMBLCID = folderChEMBL + '/PubChemCIDs.txt'  # A tab separated (tsv) file, with two elements on each line.
                                                   # The first element is the ChEMBL compound ID.
                                                   # The second element is the PubChem CID that corresponds to the ChEMBL compound ID in the first element.

    # COSMIC files used.
    cosmicData = folderCOSMIC + '/CosmicCompleteExport_v59_220512.tsv'  # The file containing all the COSMIC data.
    cosmicParsedGene = folderCOSMIC + '/COSMICParsedGene.txt'  # A file of tuples with forty-nine elements, with on tuple on each line.
                                                               # The first element is the gene ID.
                                                               # The second element is the Ensembl transcript ID that corresponds to the gene ID in the first element.
                                                               # The third element is the HGNC gene ID that corresponds to the gene ID in the first element.
                                                               # Elements four through forty-nine indicate the number of times that the gene in the first element was observed in each of the forty-six primary sites.
                                                               #     The order of the primary sites in the tuple is the same as the list primarySites in the parseCOSMIC script.
    cosmicParsedGene2Mutation = folderCOSMIC + '/COSMICParsedGene2Mutation.txt'  # A file of tuples with two elements, with one tuple on each line.
                                                                                 # The first element is the gene ID.
                                                                                 # The second element is the mutation ID.
    cosmicParsedMutation = folderCOSMIC + '/COSMICParsedMutation.txt'  # A file of tuples with sixty-five elements, with on tuple on each line.
                                                                       # The first element is the mutation ID.
                                                                       # The second element is the amino acid change.
                                                                       # The third element is 'Somatic', 'Unknown' or 'Germline' depending on whether the mutation is somatic, unknown or germline respectively.
                                                                       # Elements four through twenty indicate the type of the mutation in the first element. Only one of the sixteen options will be set to 1 (indicating the mutation is of that type), the rest will be set to 0.
                                                                       #     The order of the mutation types in the tuple is the same as the list mutationTypes in the parseCOSMIC script.
                                                                       # Elements twenty-one through sixty-five indicate the number of times that the mutation in the first element was observed in each of the forty-six primary sites.
                                                                       #     The order of the primary sites in the tuple is the same as the list primarySites in the parseCOSMIC script.

    # DrugBank files used.
    DBTargetFasta = folderDB + '/DrugBankApprovedTargets.fasta'  # A FASTA format file containing all the approved protein drug targets in DrugBank.
    DBTargetExternalLinks = folderDB + '/ExternalTargetLinks.csv'  # A csv file containing the external database cross-references for the DrugBank targets.
    DBXML = folderDB + '/DrugBank.xml'  # All the drugs in DrugBank, including their target information.
    DBDrugIDs = folderDB + '/DBDrugs.txt'  # A tab separated (tsv) file, with five elements on each line.
                                           # The first element is the DrugBank ID of the drug.
                                           # The second element is the name of the drug as recorded by DrugBank.
                                           # The third element is a semi-colon separated list of all the DrugBank drug groups that the drug is a member of.
                                           # The fourth element is the CAS number of the drug as recorded by DrugBank.
                                           # The fifth element is a semi-colon separated list of PubChem CIDs that DrugBank has linked to the drugs.
    DBTargetIDs = folderDB + '/DBTargets.txt'  # A tab separated (tsv) file, with two elements on each line.
                                               # The first element is a UniProt accession of an approved drug target.
                                               # The second element is a semi-colon separated list of all the drugs that are approved and target the protein in the first element.

    # Ensembl files used.
    ensemblGeneIDs = folderEnsembl + '/EnsemblGeneIDs.txt'  # A file with one Ensembl gene ID on each line.
    ensemblExternalIDsOne = folderEnsembl + '/ExternalOne.txt'  # Two external ID files are needed due to Ensembl Biomart limitations on the number of types of external ID that can be downloaded in one request.
    ensemblExternalIDsTwo = folderEnsembl + '/ExternalTwo.txt'  # Two external ID files are needed due to Ensembl Biomart limitations on the number of types of external ID that can be downloaded in one request.
    ensemblTranscripts = folderEnsembl + '/Transcript.txt'  # The file containing the results of extracting the transcript information from Ensembl using the Perl API.
    ensemblParsedTranscripts = folderEnsembl + '/ParsedTranscript.txt'  # A file of 6-tuples, with one on each line.
                                                                        # The first element is the Ensembl gene ID.
                                                                        # The second element is the number of different transcripts that the gene in the first element generates.
                                                                        # The third element is the number of transcripts the gene in the first element generates that are protein coding.
                                                                        # The fourth element is the number of transcripts the gene in the first element generates that are retain intron.
                                                                        # The fifth element is the number of transcripts the gene in the first element generates that are processed transcript.
                                                                        # The sixth element is the number of transcripts the gene in the first element generates that are nonsense mediated decay.
    ensemblGermSNPResults = folderEnsembl + '/GermSNPs.txt'  # The file containing the results of extracting the germline variants from Ensembl using the Perl API.
    ensemblParsedGermVariants = folderEnsembl + '/ParsedGermVariants.txt'  # A file of 25-tuples, with one on each line.
                                                                           # The first element is the Ensembl transcript ID.
                                                                           # The second element is mutation ID.
                                                                           # The third element is the Ensembl gene ID that the transcript in the first element comes from.
                                                                           # The fourth element is the amino acid change cause by the mutation ('NA' if no change is caused/known).
                                                                           # Elements five through twenty-five indicate the consequence of the mutation identified in the second element on the transcript in the first element. A 1 for an element indicates that the consequence corresponding to the element occurs.
                                                                           #     The order of the primary sites in the tuple is the same as the dictionary consequenceDict in the code.
    ensemblTaxonomyMap = folderEnsembl + '/EnsemblTaxonomyMap.txt'  # A tab separated (tsv) file, with two elements on each line.
                                                                    # The first element is the organism ID.
                                                                    # The second element is the name of the organism that corresponds to the ID in the first element.
    ensemblHomologData = folderEnsembl + '/EnsemblHomologs.txt'  # Homology data from the compara database using the standard mutli species set of species.
    ensemblGenomesTaxonomyMap = folderEnsembl + '/EnsemblGenomesTaxonomyMap.txt'  # A tab separated (tsv) file, with two elements on each line.
                                                                                  # The first element is the organism ID.
                                                                                  # The second element is the name of the organism that corresponds to the ID in the first element.
    ensemblGenomesHomologData = folderEnsembl + '/EnsemblGenomesHomologs.txt'  # Homology data from teh compara database using the pan_homology set of species.
    ensemblParsedHomology = folderEnsembl + '/EnsemblParsedHomologs.txt'  # A file containing 10-tuples, with one on each line.
                                                                          # The first element is the human Ensembl gene ID.
                                                                          # The second element is the ID of the gene that is homologus to the gene in the first element.
                                                                          # The third element is the organism that the gene in the second element comes from.
                                                                          # The fourth element is the type of homolog.
                                                                          # The fifth element is the most recent common ancestor of the genes in the first two elements.
                                                                          # The sixth element is a value for dN (non-synonymous substitutions per non-synonymous site) (-1.0 if a value is not present in Ensembl for this).
                                                                          # The seventh element is a value for dS (synonymous substitutions per synonymous site) (-1.0 if a value is not present in Ensembl for this).
                                                                          # The eighth element is the percentage of the peptide which has been aligned (-1.0 if a value is not present in Ensembl for this).
                                                                          # The ninth element is the percentage of identity between both homologs (-1.0 if a value is not present in Ensembl for this).
                                                                          # The tenth element is the percentage of positivity (similarity) between both homologs (-1.0 if a value is not present in Ensembl for this).

    # Ensembl Perl API files used.
    ensemblVariationScript = folderEnsemblPerlAPI + '/EnsemblVariant.pl'  # The perl script used to extract the variant information from Ensembl.
    ensemblHomologScript = folderEnsemblPerlAPI + '/EnsemblHomologs.pl'  # The perl script used to extract the standard multi species homology data from Ensembl Compara.
    ensemblGenomesHomologScript = folderEnsemblPerlAPI + '/EnsemblGenomesHomologs.pl'  # The perl script used to extract the pan_homology homology data from Ensembl Compara.


    # epestfind files used.
    epestfindExe = folderEpestfind + '/epestfind.exe'  # The location of the epestfind executable.

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
    pepstatsExe = folderPepstats + '/pepstats.exe'  # The location of the pepstats executable.

    # SEG files used.
    SEGExe = folderSEG + '/segmasker.exe'  # The location of the SEG executable.

    # TTD files used.
    TTDTargets = folderTTD + '/TTDTargetDataset.txt'  # The TTD drug target database in the raw format.
    TTDUPAccessions = folderTTD + '/UPAccessions.txt'  # A file of UniProt accessions, one on each line.
                                                       # Each UP accession in the file is the target of an approved drug (as recorded by the TTD).
    TTDDrugXref = folderTTD + '/TTDDrugXref.txt'  # The external cross-reference data for the drugs in the database.
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
    tableEnsemblGene = schemaProteins + '.ensemblgene'  # Table recording the information about Ensembl genes.
    tableCancerGene = schemaProteins + '.cancergene'  # A table containing information on the cancer implication status of the proteins.
    tableGOInfo = schemaProteins + '.goinfo'  # A table containing information about GO terms.
    tableProteinInfo = schemaProteins + '.proteininfo'  # A table containing the information about the proteins.
    tableBLASTResults = schemaProteins + '.blastresults'  # A table recording the results of an all-against-all BLASTing of the proteins.
    tableNonRedundant = schemaProteins + '.nonredundant'  # A table indicating the redundancy of proteins in different datasets.
    tablePPI = schemaProteins + '.ppi'  # A table containing information about the binary PPIs.
    tableGermVariants = schemaProteins + '.germvariants'  # A table containing information about germline mutations.
    tableHomologs = schemaProteins + '.homologs'  # Table recording the homologs of human proteins.
    tableUniGene = schemaProteins + '.unigene'  # Table recording the UniGene expression clusters.
    tableUniGeneTotals = schemaProteins + '.unigenetotals'  # Table recording the totals for each UniGene expression possibility.
    tableUniProt2Ensembl = schemaProteins + '.uniprot2ensembl'  # Table cross-referencing UniProt accessions with Ensembl gene, transcript and protein IDs.
    tableUniProt2GO = schemaProteins + '.uniprot2go'  # Table cross-referencing UniProt accessions and GO term IDs.
    tableUniProt2HGNC = schemaProteins + '.uniprot2hgnc'  # Table cross-referencing UniProt accessions and HGNC gene IDs.
    tableUniProt2UniGene = schemaProteins + '.uniprot2unigene'  # Table cross-referencing UniProt accessions and UniGene cluster IDs.
    tableDrugs = schemaProteins + '.drugs'  # A table containing information on the FDA approved drugs.
    tablePathways = schemaProteins + '.pathways'  # A table containing pathway information about the proteins.
    tableStability = schemaProteins + '.stability'  # A table containing stability information about the proteins.
    tableCOSMICGene = schemaProteins + '.cosmicgene'
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
    gaDatasetToGenerate = ''  # The dataset that you want to generate. Must be one of:
                              # IonChannel, GPCR, Kinase, Protease, All, CancerTarg, CancerType, CancerProt, CancerCTNCNT
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
            gaDatasetToGenerate = chunks[0]

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
            updaters.updateCancer.main(CGCParsed, cancerTargets, UPExternalLinks, UPHumanAccessionMap, TTDTarget2Drug, DBTargetIDs,
                                       schemaProteins, tableCancerGene, DATABASEPASSWORD)
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
                                                    ChEMBLUPAccessions, UPHumanAccessionMap, UPDrugIDs, folderCulling,
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
        # Generate the output files for the dataset generation.
        # Assumes that the columns are identical for both tables/views.
        print '\nGenerating GA Data File.'
        outputDirectory = folderGAData + '/' + gaDatasetToGenerate
        if os.path.isdir(outputDirectory):
            shutil.rmtree(outputDirectory)
        os.mkdir(outputDirectory)
        nonRedundantGOOutputDirectory = outputDirectory + '/GO-NonRedundant'
        if os.path.isdir(nonRedundantGOOutputDirectory):
            shutil.rmtree(nonRedundantGOOutputDirectory)
        os.mkdir(nonRedundantGOOutputDirectory)
        datasetFastaFile = outputDirectory + '/' + gaDatasetToGenerate + '.fasta'
        gaNonRedundantData = outputDirectory + '/' + gaDatasetToGenerate + '-NonRedundant.txt'
        columnNRDataLocation = outputDirectory + '/Columns-NonRedundant.txt'
        ECNRDataLocation = outputDirectory + '/ECNumbers-NonRedundant.txt'
        subcellNRLocation = outputDirectory + '/SubcellularLocation-NonRedundant.txt'
        healthStateNRLocation = outputDirectory + '/HealthState-NonRedundant.txt'
        bodySiteNRLocation = outputDirectory + '/BodySite-NonRedundant.txt'
        developmentalStageNRLocation = outputDirectory + '/DevelopmentalStage-NonRedundant.txt'
        gaRedundantData = outputDirectory + '/' + gaDatasetToGenerate + '-Redundant.txt'
        columnRDataLocation = outputDirectory + '/Columns-Redundant.txt'
        ECRDataLocation = outputDirectory + '/ECNumbers-Redundant.txt'
        subcellRLocation = outputDirectory + '/SubcellularLocation-Redundant.txt'
        healthStateRLocation = outputDirectory + '/HealthState-Redundant.txt'
        bodySiteRLocation = outputDirectory + '/BodySite-Redundant.txt'
        developmentalStageRLocation = outputDirectory + '/DevelopmentalStage-Redundant.txt'
        blastOutput = outputDirectory + '/BlastData.txt'
        gaAllData = outputDirectory + '/' + gaDatasetToGenerate + '-All.txt'
        
        positiveNRViewDict = {'All' : viewAllAllTargNRP,
                            'GPCR' : viewTypeGPCRTargNRP,
                            'IonChannel' : viewTypeIonTargNRP,
                            'Kinase' : viewTypeKinaseTargNRP,
                            'Protease' : viewTypeProteaseTargNRP,
                            'CancerTarg' : viewIllCancerTargNRP,
                            'CancerType' : viewIllCancerTypeNRP,
                            'CancerProt' : viewIllCancerProtNRP,
                            'CancerCTNCNT' : viewIllCancerCTNCNTNRP}
        unlabelledNRViewDict = {'All' : viewAllAllTargNRN,
                            'GPCR' : viewTypeGPCRTargNRN,
                            'IonChannel' : viewTypeIonTargNRN,
                            'Kinase' : viewTypeKinaseTargNRN,
                            'Protease' : viewTypeProteaseTargNRN,
                            'CancerTarg' : viewIllCancerTargNRN,
                            'CancerType' : viewIllCancerTypeNRN,
                            'CancerProt' : viewIllCancerProtNRN,
                            'CancerCTNCNT' : viewIllCancerCTNCNTNRN}
        positiveRedundantViewSQLDict = {'All' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'Y\'',
                            'GPCR' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'Y\' AND ModeOfAction = \'G-protein coupled receptor\'',
                            'IonChannel' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'Y\' AND ModeOfAction = \'Ion Channel\'',
                            'Kinase' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'Y\' AND ModeOfAction = \'Kinase\'',
                            'Protease' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'Y\' AND ModeOfAction = \'Protease\'',
                            'CancerTarg' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                            'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                            'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'Y\' AND '
                                                'cancer.Target=\'Y\''
                                            ),
                            'CancerType' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                            'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                            'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'Y\' AND '
                                                'cancer.Target=\'Y\''
                                            ),
                            'CancerProt' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                            'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                            'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'Y\''
                                            ),
                            'CancerCTNCNT' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                            'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                            'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'Y\' AND '
                                                'cancer.Target=\'Y\''
                                            ),}
        unlabelledRedundantViewSQLDict = {'All' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'N\'',
                            'GPCR' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'N\' AND ModeOfAction = \'G-protein coupled receptor\'',
                            'IonChannel' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'N\' AND ModeOfAction = \'Ion Channel\'',
                            'Kinase' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'N\' AND ModeOfAction = \'Kinase\'',
                            'Protease' : 'SELECT UPAccession, Sequence FROM ' + tableProteinInfo + ' WHERE Target=\'N\' AND ModeOfAction = \'Protease\'',
                            'CancerTarg' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                            'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                            'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'Y\' AND '
                                                'cancer.Target=\'N\''
                                            ),
                            'CancerType' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                             'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                             'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'N\' AND '
                                                'prot.Target=\'Y\''
                                            ),
                            'CancerProt' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                             'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                             'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'N\''
                                            ),
                            'CancerCTNCNT' : ('SELECT '
                                                'prot.UPAccession, prot.Sequence '
                                             'FROM ' +
                                                tableProteinInfo + ' AS prot, ' +
                                                tableCancerGene + ' AS cancer '
                                             'WHERE '
                                                'cancer.UPAccession = prot.UPAccession AND '
                                                'cancer.Cancer=\'N\' AND '
                                                'prot.Target=\'N\''
                                            )}
        positiveColumnNameDict = {'All' : 'AllTargetPositive',
                            'GPCR' : 'GPCRTargetPositive',
                            'IonChannel' : 'IonChannelTargetPositive',
                            'Kinase' : 'KinaseTargetPositive',
                            'Protease' : 'ProteaseTargetPositive',
                            'CancerTarg' : 'CancerTargetPositive',
                            'CancerType' : 'CancerTypePositive',
                            'CancerProt' : 'CancerProteinPositive',
                            'CancerCTNCNT' : 'CancerCTNCNTPositive'}
        unlabelledColumnNameDict = {'All' : 'AllTargetNegative',
                            'GPCR' : 'GPCRTargetNegative',
                            'IonChannel' : 'IonChannelTargetNegative',
                            'Kinase' : 'KinaseTargetNegative',
                            'Protease' : 'ProteaseTargetNegative',
                            'CancerTarg' : 'CancerTargetNegative',
                            'CancerType' : 'CancerTypeNegative',
                            'CancerProt' : 'CancerProteinNegative',
                            'CancerCTNCNT' : 'CancerCTNCNTNegative'}
        positiveColumn = positiveColumnNameDict[gaDatasetToGenerate]  # Get the name of the positive observation column in tableNonRedundant.
        unlabelledColumn = unlabelledColumnNameDict[gaDatasetToGenerate]  # Get the name of the unlabelled observation column in tableNonRedundant.
        positiveNRView = positiveNRViewDict[gaDatasetToGenerate]  # Get the name of the positive non-redundant observation view.
        unlabelledNRView = unlabelledNRViewDict[gaDatasetToGenerate]  # Get the name of the unlabelled non-redundant observation view.
        positiveRedundantViewSQL = positiveRedundantViewSQLDict[gaDatasetToGenerate]  # Get the SQL for generating the positive redundant observation view.
        unlabelledRedundantViewSQL = unlabelledRedundantViewSQLDict[gaDatasetToGenerate]  # Get the SQL for generating the unlabelled redundant observation view.

        # Generate a FASTA format file of all the proteins in the dataset.
        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        positiveProteins = cursor.execute(positiveRedundantViewSQL)
        positiveProteins = cursor.fetchall()
        unlabelledProteins = cursor.execute(unlabelledRedundantViewSQL)
        unlabelledProteins = cursor.fetchall()
        results = positiveProteins + unlabelledProteins
        utilities.list2file.main(['>' + '\n'.join([j[0], j[1]]) for j in results], datasetFastaFile)
        mysql.closeConnection(conn, cursor)

        # Ideally you would just SELECT * FROM viewToGenerateFrom. However, the view is so complex to generate that it is substantially quicker to
        # extract the data from the the tables themselves, and then generate the view's data manually. (At least for me it is)
        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        # Get the names of all the columns in the view of the non-redundant dataset, and therefore the names of the features in the dataset.
        cursor.execute('SHOW COLUMNS FROM ' + positiveNRView)
        columns = cursor.fetchall()
        columns = [i[0] for i in columns]

        # Get the name of all the columns in the non-redundant dataset table.
        nonredColumns = cursor.execute('SHOW COLUMNS FROM ' + tableNonRedundant)
        nonredColumns = cursor.fetchall()
        nonredColumns = [i[0] for i in nonredColumns]
        # Get the data from the non-redundant dataset table.
        nonredundantPositiveProteins = cursor.execute('SELECT UPAccession FROM ' + tableNonRedundant + ' WHERE ' + positiveColumn + '="Y"')
        nonredundantPositiveProteins = cursor.fetchall()
        positiveNonRedundantProteinAccs = [i[0] for i in nonredundantPositiveProteins]  # Determine the UniProt accessions of the proteins in the positive non-redundant dataset.
        nonredundantUnlabelledProteins = cursor.execute('SELECT UPAccession FROM ' + tableNonRedundant + ' WHERE ' + unlabelledColumn + '="Y"')
        nonredundantUnlabelledProteins = cursor.fetchall()
        unlabelledNonRedundantProteinAccs = [i[0] for i in nonredundantUnlabelledProteins]  # Determine the UniProt accessions of the proteins in the unlabelled non-redundant dataset.
        redundantPositiveProteins = cursor.execute(positiveRedundantViewSQL)
        redundantPositiveProteins = cursor.fetchall()
        positiveRedundantProteinAccs = [i[0] for i in redundantPositiveProteins if i[0] not in positiveNonRedundantProteinAccs]  # Determine the UniProt accessions of the proteins in the positive non-redundant dataset.
        redundantUnlabelledProteins = cursor.execute(unlabelledRedundantViewSQL)
        redundantUnlabelledProteins = cursor.fetchall()
        unlabelledRedundantProteinAccs = [i[0] for i in redundantUnlabelledProteins if i[0] not in unlabelledNonRedundantProteinAccs]  # Determine the UniProt accessions of the proteins in the unlabelled non-redundant dataset.

        # Get the names of the columns in the UniProt protein information table.
        protColumns = cursor.execute('SHOW COLUMNS FROM ' + tableProteinInfo)
        protColumns = cursor.fetchall()
        protColumns = [i[0] for i in protColumns]
        # Get the data in the UniProt protien information table.
        proteins = cursor.execute('SELECT * FROM ' + tableProteinInfo)
        proteins = cursor.fetchall()
        positiveNonRedundantProteins = dict([(i[0], i) for i in proteins if i[0] in positiveNonRedundantProteinAccs])  # Get the UniProt protein information for the proteins in the non-redundant positive dataset.
        unlabelledNonRedundantProteins = dict([(i[0], i) for i in proteins if i[0] in unlabelledNonRedundantProteinAccs])  # Get the UniProt protein information for the proteins in the non-redundant unlabelled dataset.
        positiveRedundantProteins = dict([(i[0], i) for i in proteins if i[0] in positiveRedundantProteinAccs])  # Get the UniProt protein information for the proteins in the redundant positive dataset.
        unlabelledRedundantProteins = dict([(i[0], i) for i in proteins if i[0] in unlabelledRedundantProteinAccs])  # Get the UniProt protein information for the proteins in the redundant unlabelled dataset.

        # Get the names of the columns from the stability information table.
        stabilityColumns = cursor.execute('SHOW COLUMNS FROM ' + tableStability)
        stabilityColumns = cursor.fetchall()
        stabilityColumns = [i[0] for i in stabilityColumns]
        # Get the data from teh stability information table.
        stability = cursor.execute('SELECT * FROM ' + tableStability)
        stability = cursor.fetchall()
        positiveNonRedundantStability = dict([(i[0], i) for i in stability if i[0] in positiveNonRedundantProteinAccs])  # Get the stability information for the proteins in the positive dataset.
        unlabelledNonRedundantStability = dict([(i[0], i) for i in stability if i[0] in unlabelledNonRedundantProteinAccs])  # Get the stability information for the proteins int he unlabelled dataset.
        positiveRedundantStability = dict([(i[0], i) for i in stability if i[0] in positiveRedundantProteinAccs])  # Get the stability information for the proteins in the positive dataset.
        unlabelledRedundantStability = dict([(i[0], i) for i in stability if i[0] in unlabelledRedundantProteinAccs])  # Get the stability information for the proteins int he unlabelled dataset.

        # Get the data from the PPI information views.
        ppi = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_ppi')
        ppi = cursor.fetchall()
        positiveNonRedundantPPI = dict([(i[0], i) for i in ppi if i[0] in positiveNonRedundantProteinAccs])  # Get the PPI information for the proteins in the non-redundant positive dataset.
        unlabelledNonRedundantPPI = dict([(i[0], i) for i in ppi if i[0] in unlabelledNonRedundantProteinAccs])  # Get the PPI information for the proteins in the non-redundant positive dataset.
        positiveRedundantPPI = dict([(i[0], i) for i in ppi if i[0] in positiveRedundantProteinAccs])  # Get the PPI information for the proteins in the redundant positive dataset.
        unlabelledRedundantPPI = dict([(i[0], i) for i in ppi if i[0] in unlabelledRedundantProteinAccs])  # Get the PPI information for the proteins in the redundant positive dataset.

        # Get the names of the columns in the expression information view.
        expressionColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_expression')
        expressionColumns = cursor.fetchall()
        expressionColumns = [i[0] for i in expressionColumns]
        # Get the data from teh expression information view.
        expression = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_expression')
        expression = cursor.fetchall()
        positiveNonRedundantExpression = dict([(i[0], i) for i in expression if i[0] in positiveNonRedundantProteinAccs])  # Get the expression information for the proteins in the non-redundant positive dataset.
        unlabelledNonRedundantExpression = dict([(i[0], i) for i in expression if i[0] in unlabelledNonRedundantProteinAccs])  # Get the expression information for the proteins in the non-redundant positive dataset.
        positiveRedundantExpression = dict([(i[0], i) for i in expression if i[0] in positiveRedundantProteinAccs])  # Get the expression information for the proteins in the redundant positive dataset.
        unlabelledRedundantExpression = dict([(i[0], i) for i in expression if i[0] in unlabelledRedundantProteinAccs])  # Get the expression information for the proteins in the redundant positive dataset.

        # Get the names of the columns from the paralogs view.
        paralogColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_paralogs')
        paralogColumns = cursor.fetchall()
        paralogColumns = [i[0] for i in paralogColumns]
        # Get the data from the paralogs view.
        paralog = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_paralogs')
        paralog = cursor.fetchall()
        positiveNonRedundantParalog = dict([(i[0], i) for i in paralog if i[0] in positiveNonRedundantProteinAccs])  # Get the paralog information for the proteins in the non-redundant positive dataset.
        unlabelledNonRedundantParalog = dict([(i[0], i) for i in paralog if i[0] in unlabelledNonRedundantProteinAccs])  # Get the paralog information for the proteins in the non-redundant positive dataset.
        positiveRedundantParalog = dict([(i[0], i) for i in paralog if i[0] in positiveRedundantProteinAccs])  # Get the paralog information for the proteins in the redundant positive dataset.
        unlabelledRedundantParalog = dict([(i[0], i) for i in paralog if i[0] in unlabelledRedundantProteinAccs])  # Get the paralog information for the proteins in the redundant positive dataset.

        # Get the column names from the transcript view.
        transcriptColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_transcripts')
        transcriptColumns = cursor.fetchall()
        transcriptColumns = [i[0] for i in transcriptColumns]
        # Get the data from the transcript view.
        transcript = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_transcripts')
        transcript = cursor.fetchall()
        positiveNonRedundantTranscript = dict([(i[0], i) for i in transcript if i[0] in positiveNonRedundantProteinAccs])  # Get the transcript information for the proteins in the non-redundant positive dataset.
        unlabelledNonRedundantTranscript = dict([(i[0], i) for i in transcript if i[0] in unlabelledNonRedundantProteinAccs])  # Get the transcript information for the proteins in the non-redundant positive dataset.
        positiveRedundantTranscript = dict([(i[0], i) for i in transcript if i[0] in positiveRedundantProteinAccs])  # Get the transcript information for the proteins in the redundant positive dataset.
        unlabelledRedundantTranscript = dict([(i[0], i) for i in transcript if i[0] in unlabelledRedundantProteinAccs])  # Get the transcript information for the proteins in the redundant positive dataset.

        # Get the column names from the germline mutations view.
        variantColumns = cursor.execute('SHOW COLUMNS FROM ' + schemaProteins + '.upacc_germvariants')
        variantColumns = cursor.fetchall()
        variantColumns = [i[0] for i in paralogColumns]
        # Get the data from the germline mutations view.
        variant = cursor.execute('SELECT * FROM ' + schemaProteins + '.upacc_germvariants')
        variant = cursor.fetchall()
        positiveNonRedundantVariant = {}
        unlabelledNonRedundantVariant = {}
        positiveRedundantVariant = {}
        unlabelledRedundantVariant = {}
        for i in variant:
            # For every mutation.
            if i[0] in positiveNonRedundantProteinAccs:
                # Get the germline mutation information for the proteins in the non-redundant positive dataset.
                if positiveNonRedundantVariant.has_key(i[0]):
                    positiveNonRedundantVariant[i[0]]['3untrans'] += i[2]
                    positiveNonRedundantVariant[i[0]]['5untrans'] += i[3]
                    positiveNonRedundantVariant[i[0]]['nonsynon'] += i[4]
                    positiveNonRedundantVariant[i[0]]['synon'] += i[5]
                else:
                    positiveNonRedundantVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}
            elif i[0] in unlabelledNonRedundantProteinAccs:
                # Get the germline mutation information for the proteins in the non-redundant positive dataset.
                if unlabelledNonRedundantVariant.has_key(i[0]):
                    unlabelledNonRedundantVariant[i[0]]['3untrans'] += i[2]
                    unlabelledNonRedundantVariant[i[0]]['5untrans'] += i[3]
                    unlabelledNonRedundantVariant[i[0]]['nonsynon'] += i[4]
                    unlabelledNonRedundantVariant[i[0]]['synon'] += i[5]
                else:
                    unlabelledNonRedundantVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}

            if i[0] in positiveRedundantProteinAccs:
                # Get the germline mutation information for the proteins in the non-redundant positive dataset.
                if positiveRedundantVariant.has_key(i[0]):
                    positiveRedundantVariant[i[0]]['3untrans'] += i[2]
                    positiveRedundantVariant[i[0]]['5untrans'] += i[3]
                    positiveRedundantVariant[i[0]]['nonsynon'] += i[4]
                    positiveRedundantVariant[i[0]]['synon'] += i[5]
                else:
                    positiveRedundantVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}
            elif i[0] in unlabelledRedundantProteinAccs:
                # Get the germline mutation information for the proteins in the non-redundant positive dataset.
                if unlabelledRedundantVariant.has_key(i[0]):
                    unlabelledRedundantVariant[i[0]]['3untrans'] += i[2]
                    unlabelledRedundantVariant[i[0]]['5untrans'] += i[3]
                    unlabelledRedundantVariant[i[0]]['nonsynon'] += i[4]
                    unlabelledRedundantVariant[i[0]]['synon'] += i[5]
                else:
                    unlabelledRedundantVariant[i[0]] = {'3untrans' : i[2], '5untrans' : i[3], 'nonsynon' : i[4], 'synon' : i[5]}
        mysql.closeConnection(conn, cursor)

        resultsNonRedundantPositive = []  # Holds the data tuples for the non-redundant positive observations.
        resultsNonRedundantUnlabelled = []  # Holds the data tuples for the non-redudnant unlabelled observations.
        resultsRedundantPositive = []  # Holds the data tuples for the redundant positive observations.
        resultsRedundantUnlabelled = []  # Holds the data tuples for the redudnant unlabelled observations.
        for i in [[positiveNonRedundantProteinAccs, 'NRP'], [positiveRedundantProteinAccs, 'RP'], [unlabelledNonRedundantProteinAccs, 'NRU'],
                  [unlabelledRedundantProteinAccs, 'RU']]:
            for j in i[0]:
                currentRecord = []  # Hold the current protein's data tuple.

                # Put the protein properties into the current protein's data tuple.
                for k in ['UPAccession', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                              'NegativelyCharged', 'PositivelyCharged', 'Basic', 'Charged', 'Polar', 'NonPolar', 'Aromatic', 'Aliphatic', 'Small', 'Tiny',
                              'PESTMotif', 'LowComplexity', 'Hydrophobicity', 'Isoelectric', 'ECNumber', 'OGlycosylation', 'NGlycosylation', 'Phosphoserine',
                              'Phosphothreonine', 'Phosphotyrosine', 'SubcellularLocation', 'TopologicalDomain', 'PredictedSubcellularLocation',
                              'SignalPeptide', 'TransmembraneHelices', 'Turns', 'AlphaHelices', 'BetaStrands', 'PredictedAlphaHelices', 'PredictedBetaSheets', 'Sequence']:
                    if i[1] == 'NRP':
                        currentRecord.append(positiveNonRedundantProteins[j][protColumns.index(k)])
                    elif i[1] == 'NRU':
                        currentRecord.append(unlabelledNonRedundantProteins[j][protColumns.index(k)])
                    elif i[1] == 'RP':
                        currentRecord.append(positiveRedundantProteins[j][protColumns.index(k)])
                    elif i[1] == 'RU':
                        currentRecord.append(unlabelledRedundantProteins[j][protColumns.index(k)])

                # Put the expression information into the current protein's data tuple.
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
                    if i[1] == 'NRP':
                        currentRecord.append(int(positiveNonRedundantExpression[j][expressionColumns.index(k)]))
                    elif i[1] == 'NRU':
                        currentRecord.append(int(unlabelledNonRedundantExpression[j][expressionColumns.index(k)]))
                    elif i[1] == 'RP':
                        currentRecord.append(int(positiveRedundantExpression[j][expressionColumns.index(k)]))
                    elif i[1] == 'RU':
                        currentRecord.append(int(unlabelledRedundantExpression[j][expressionColumns.index(k)]))

                # Put the germline mutation information into the current protein's data tuple.
                if i[1] == 'NRP':
                    if positiveNonRedundantVariant.has_key(j):
                        currentRecord.append(positiveNonRedundantVariant[j]['3untrans'])
                        currentRecord.append(positiveNonRedundantVariant[j]['5untrans'])
                        currentRecord.append(positiveNonRedundantVariant[j]['nonsynon'])
                        currentRecord.append(positiveNonRedundantVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                elif i[1] == 'NRU':
                    if unlabelledNonRedundantVariant.has_key(j):
                        currentRecord.append(unlabelledNonRedundantVariant[j]['3untrans'])
                        currentRecord.append(unlabelledNonRedundantVariant[j]['5untrans'])
                        currentRecord.append(unlabelledNonRedundantVariant[j]['nonsynon'])
                        currentRecord.append(unlabelledNonRedundantVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                elif i[1] == 'RP':
                    if positiveRedundantVariant.has_key(j):
                        currentRecord.append(positiveRedundantVariant[j]['3untrans'])
                        currentRecord.append(positiveRedundantVariant[j]['5untrans'])
                        currentRecord.append(positiveRedundantVariant[j]['nonsynon'])
                        currentRecord.append(positiveRedundantVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                elif i[1] == 'RU':
                    if unlabelledRedundantVariant.has_key(j):
                        currentRecord.append(unlabelledRedundantVariant[j]['3untrans'])
                        currentRecord.append(unlabelledRedundantVariant[j]['5untrans'])
                        currentRecord.append(unlabelledRedundantVariant[j]['nonsynon'])
                        currentRecord.append(unlabelledRedundantVariant[j]['synon'])
                    else:
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)
                        currentRecord.append(0)

                # Put the paralog information into the current protein's data tuple.
                if i[1] == 'NRP':
                    if positiveNonRedundantParalog.has_key(j):
                        currentRecord.append(int(positiveNonRedundantParalog[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'NRU':
                    if unlabelledNonRedundantParalog.has_key(j):
                        currentRecord.append(int(unlabelledNonRedundantParalog[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'RP':
                    if positiveRedundantParalog.has_key(j):
                        currentRecord.append(int(positiveRedundantParalog[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'RU':
                    if unlabelledRedundantParalog.has_key(j):
                        currentRecord.append(int(unlabelledRedundantParalog[j][-1]))
                    else:
                        currentRecord.append(0)

                # Put the PPI information into the current protein's data tuple.
                if i[1] == 'NRP':
                    if positiveNonRedundantPPI.has_key(j):
                        currentRecord.append(int(positiveNonRedundantPPI[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'NRU':
                    if unlabelledNonRedundantPPI.has_key(j):
                        currentRecord.append(int(unlabelledNonRedundantPPI[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'RP':
                    if positiveRedundantPPI.has_key(j):
                        currentRecord.append(int(positiveRedundantPPI[j][-1]))
                    else:
                        currentRecord.append(0)
                elif i[1] == 'RU':
                    if unlabelledRedundantPPI.has_key(j):
                        currentRecord.append(int(unlabelledRedundantPPI[j][-1]))
                    else:
                        currentRecord.append(0)

                # Put the transcript information into the current protein's data tuple.
                if i[1] == 'NRP':
                    if positiveNonRedundantTranscript.has_key(j):
                        currentRecord.append(int(positiveNonRedundantTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)
                elif i[1] == 'NRU':
                    if unlabelledNonRedundantTranscript.has_key(j):
                        currentRecord.append(int(unlabelledNonRedundantTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)
                elif i[1] == 'RP':
                    if positiveRedundantTranscript.has_key(j):
                        currentRecord.append(int(positiveRedundantTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)
                elif i[1] == 'RU':
                    if unlabelledRedundantTranscript.has_key(j):
                        currentRecord.append(int(unlabelledRedundantTranscript[j][transcriptColumns.index('ProteinCodingTranscripts')]))
                    else:
                        currentRecord.append(1)

                # Put the stability information into the current protein's data tuple.
                if i[1] == 'NRP':
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(positiveNonRedundantStability[j][stabilityColumns.index(k)])
                elif i[1] == 'NRU':
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(unlabelledNonRedundantStability[j][stabilityColumns.index(k)])
                elif i[1] == 'RP':
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(positiveRedundantStability[j][stabilityColumns.index(k)])
                elif i[1] == 'RU':
                    for k in ['HalfLife', 'InstabilityIndex']:
                        currentRecord.append(unlabelledRedundantStability[j][stabilityColumns.index(k)])

                if i[1] == 'NRP':
                    resultsNonRedundantPositive.append(tuple(currentRecord))
                elif i[1] == 'NRU':
                    resultsNonRedundantUnlabelled.append(tuple(currentRecord))
                elif i[1] == 'RP':
                    resultsRedundantPositive.append(tuple(currentRecord))
                elif i[1] == 'RU':
                    resultsRedundantUnlabelled.append(tuple(currentRecord))

        # Generate the dataset and any additional information about the proteins in it (expression information, EC numbers, etc.).
        utilities.gadatageneration.pulearning(resultsNonRedundantPositive, resultsNonRedundantUnlabelled, columns, gaNonRedundantData,
                                              columnNRDataLocation, ECNRDataLocation, subcellNRLocation,
                                              healthStateNRLocation, bodySiteNRLocation, developmentalStageNRLocation)
        utilities.gadatageneration.pulearning(resultsRedundantPositive, resultsRedundantUnlabelled, columns, gaRedundantData,
                                              columnRDataLocation, ECRDataLocation, subcellRLocation,
                                              healthStateRLocation, bodySiteRLocation, developmentalStageRLocation)

        # Generate the entire dataset.
        nonRedundantData = set([])
        readNonRedundant = open(gaNonRedundantData, 'r')
        headerOne = readNonRedundant.readline()
        for line in readNonRedundant:
            nonRedundantData.add(line)
        readNonRedundant.close()
        redundantData = set([])
        readRedundant = open(gaRedundantData, 'r')
        headerOne = readRedundant.readline()
        for line in readRedundant:
            redundantData.add(line)
        readRedundant.close()
        allData = redundantData | nonRedundantData
        writeAll = open(gaAllData, 'w')
        writeAll.write(headerOne)
        for i in allData:
            writeAll.write(i)
        writeAll.close()

        # Determine the BLAST similarity info for the proteins in the non-redundant dataset.
        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        allRedundantProteinAccs = positiveRedundantProteinAccs + unlabelledRedundantProteinAccs
        blastTableQuery = cursor.execute('SELECT * FROM ' + tableBLASTResults + ' WHERE ProteinA IN (' + ','.join(['\'' + i + '\'' for i in allRedundantProteinAccs]) + ') AND ProteinB IN (' + ','.join(['\'' + i + '\'' for i in allRedundantProteinAccs]) + ')')
        blastTableQuery = cursor.fetchall()
        blastSimilarities = dict([(i, {'Protein' : [], 'Similarity' : []}) for i in allRedundantProteinAccs])
        for i in blastTableQuery:
            protA = i[0]
            protB = i[1]
            similarity = i[2]
            if similarity <= 20:
                # If the similarity of the proteins is below the 20% threshold, then ignore the similarity.
                continue
            blastSimilarities[protB]['Protein'].append(protA)
            blastSimilarities[protB]['Similarity'].append(similarity)
            blastSimilarities[protA]['Protein'].append(protB)
            blastSimilarities[protA]['Similarity'].append(similarity)
        mysql.closeConnection(conn, cursor)
        writeOut = open(blastOutput, 'w')
        for i in blastSimilarities:
            if blastSimilarities[i]['Protein']:
                similarities, proteins = zip(*sorted(zip(blastSimilarities[i]['Similarity'], blastSimilarities[i]['Protein']), reverse=True))
            else:
                similarities = []
                proteins = []
            writeOut.write(i)
            writeOut.write('\t')
            writeOut.write(','.join(proteins))
            writeOut.write('\t')
            writeOut.write(','.join([str(i) for i in similarities]))
            writeOut.write('\n')
        writeOut.close()

        # Determine the GO term information for the proteins in the dataset.
        conn, cursor = mysql.openConnection(DATABASEPASSWORD, schemaProteins)
        # Get the GO term data for the positive proteins.
        cursor = mysql.tableSELECT(cursor, '*', tableUniProt2GO, 'UPAccession IN ("' + '","'.join(positiveNonRedundantProteinAccs) + '")')
        results = cursor.fetchall()
        positiveNonRedundantGOIDs = {}
        for i in results:
            UPAccession = i[0]
            GOTermID = i[1]
            if positiveNonRedundantGOIDs.has_key(UPAccession):
                positiveNonRedundantGOIDs[UPAccession].add(GOTermID)
            else:
                positiveNonRedundantGOIDs[UPAccession] = set([GOTermID])
        positiveNonRedundantGOTerms = {}
        for i in positiveNonRedundantGOIDs.keys():
            UPAccession = i
            positiveNonRedundantGOTerms[UPAccession] = {'biological_process' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                          'cellular_component' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                          'molecular_function' : {'LevelOne' : set([]), 'LevelTwo' : set([])}
                                          }
            GOIDs = [str(i) for i in positiveNonRedundantGOIDs[i]]
            cursor = mysql.tableSELECT(cursor, 'GOType, LevelOne, LevelTwo', tableGOInfo, 'GOTermID IN ("' + '","'.join(GOIDs) + '")')
            results = cursor.fetchall()
            for j in results:
                GOType = j[0]
                levelOne = j[1].split(';')
                levelTwo = j[2].split(';')

                if levelOne != ['NA']:
                    positiveNonRedundantGOTerms[UPAccession][GOType]['LevelOne'].update(levelOne)
                if levelTwo != ['NA']:
                    positiveNonRedundantGOTerms[UPAccession][GOType]['LevelTwo'].update(levelTwo)

        # Get the GO term data for the unlabelled proteins.
        cursor = mysql.tableSELECT(cursor, '*', tableUniProt2GO, 'UPAccession IN ("' + '","'.join(unlabelledNonRedundantProteinAccs) + '")')
        results = cursor.fetchall()
        unlabelledNonRedundantGOID = {}
        for i in results:
            UPAccession = i[0]
            GOTermID = i[1]
            if unlabelledNonRedundantGOID.has_key(UPAccession):
                unlabelledNonRedundantGOID[UPAccession].add(GOTermID)
            else:
                unlabelledNonRedundantGOID[UPAccession] = set([GOTermID])
        unlabelledNonRedundantGOTerms = {}
        for i in unlabelledNonRedundantGOID.keys():
            UPAccession = i
            unlabelledNonRedundantGOTerms[UPAccession] = {'biological_process' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                             'cellular_component' : {'LevelOne' : set([]), 'LevelTwo' : set([])},
                                             'molecular_function' : {'LevelOne' : set([]), 'LevelTwo' : set([])}
                                             }
            GOIDs = [str(i) for i in unlabelledNonRedundantGOID[i]]
            cursor = mysql.tableSELECT(cursor, 'GOType, LevelOne, LevelTwo', tableGOInfo, 'GOTermID IN ("' + '","'.join(GOIDs) + '")')
            results = cursor.fetchall()
            for j in results:
                GOType = j[0]
                levelOne = j[1].split(';')
                levelTwo = j[2].split(';')

                if levelOne != ['NA']:
                    unlabelledNonRedundantGOTerms[UPAccession][GOType]['LevelOne'].update(levelOne)
                if levelTwo != ['NA']:
                    unlabelledNonRedundantGOTerms[UPAccession][GOType]['LevelTwo'].update(levelTwo)
        mysql.closeConnection(conn, cursor)

        # Generate the summary GO term information for the positive and negative observations in the dataset.
        utilities.generateGOsummarydata.main(positiveNonRedundantGOTerms, unlabelledNonRedundantGOTerms, nonRedundantGOOutputDirectory)

if __name__ == '__main__':
    main(sys.argv[1:])