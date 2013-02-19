'''
Created on 11 Oct 2011

@author: Simon Bull
'''

import time
import urllib2
import re

import utilities.MySQLaccess as mysql
import utilities.list2file

def main(ChEMBLTargets, ChEMBLCID, databasePassword, ChEMBLSchema):
    """
    Returns two files.
        ChEMBLTargets contains mapping of approved drug information to the UniProt accessions of the proteins that the drug targets.
        ChEMBLCID contains mappings of ChEMBL compound IDs to PubChem CIDs.
    ChEMBLTargets - A tab separated (tsv) file, with seven elements on each line. There is one line for each protein that an approved target targets.
        The first element is the UniProt accession for the target protein.
        The second element is ChEMBL ID of the compound.
        The third element is the name of the protein in the first element.
        The fourth element is the activity relation (=, <, <=, > or >=).
        The fifth element is the value of the activity observed between the compound and the protein.
        The sixth element is the units of the value in the fifth element.
        The seventh element is the type of activity measured (Ki, Kd, EC50, etc.)
    ChEMBLCID - A tab separated (tsv) file, with two elements on each line.
        The first element is the ChEMBL compound ID.
        The second element is the PubChem CID that corresponds to the ChEMBL compound ID in the first element.
    """

    # Create the query for extracting targets of approved drugs.
    queryTarget2Compound = """
    SELECT
        td.protein_accession,
        md.molregno,
        md.chembl_id,
        td.pref_name,
        act.relation,
        act.standard_value,
        act.standard_units,
        act.standard_type
    FROM
        molecule_dictionary md,
        activities act,
        assays a,
        assay2target a2t,
        target_dictionary td
    WHERE
        md.max_phase = '4'
        and md.molecule_type = "Small molecule"
        and md.molregno = act.molregno
        and act.assay_id = a.assay_id
        and a.assay_id = a2t.assay_id
        and a2t.confidence_score >= '4'
        and a2t.tid = td.tid
        and td.target_type = 'PROTEIN'
        and td.db_source = 'SWISS-PROT'
    """

    # Connect to the ChEMBL schema in the database, and run the target/compound extraction query.
    conn, cursor = mysql.openConnection(databasePassword, ChEMBLSchema)
    cursor.execute(queryTarget2Compound)
    resultsTarget2Compound = cursor.fetchall()
    resultsTarget2Compound = list(set(resultsTarget2Compound))

    # Generate the output file of UniProt accession to drug information mappings.
    target2CompoundsDict = {}
    chemblID2Molregno = {}
    ChEMBLIDs = set([])
    writeOut = open(ChEMBLTargets, 'w')
    for i in resultsTarget2Compound:
        chemblID2Molregno[str(i[2])] = str(i[1])
        ChEMBLIDs.add(str(i[2]))
        writeOut.write(i[0] + '\t' + str(i[1]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i[5]) + '\t' + str(i[6]) + '\t' + str(i[7]) + '\n')
    writeOut.close()
    mysql.closeConnection(conn, cursor)

    # Use Entrez EUtils to get the PubChem CIDs for the ChEMBL IDs.
    molregno2CIDs = set([])
    for i in ChEMBLIDs:
        fetchURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term=' + i
        req = urllib2.Request(fetchURL)
        fetchResult = urllib2.urlopen(req)
        output = fetchResult.read()
        result = re.search('(?<=<Id>)[0-9]*(?=</Id>)', output)
        if result:
            result = result.group(0)
            molregno2CIDs.add(chemblID2Molregno[i] + '\t' + result)
        fetchResult.close()
        # Sleep for 0.6 seconds in order to stay within NCBI's guidelines of not submitting more than 3 HTTP requests a second
        time.sleep(0.6)

    utilities.list2file.main(list(molregno2CIDs), ChEMBLCID)