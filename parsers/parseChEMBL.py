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
    """Extracts the UniProt IDs stored in the ChEMBL database.

    """

    extract_targets(ChEMBLTargets, ChEMBLCID, databasePassword, ChEMBLSchema)    
##    molecules = extract_targets(ChEMBLTargets, ChEMBLCID, databasePassword, ChEMBLSchema)
##    activities = extract_drugs(molecules, ChEMBLDrugs, databasePassword, ChEMBLSchema)
##    extract_activities(activities, ChEMBLActivities, databasePassword, ChEMBLSchema)


def extract_targets(ChEMBLTargets, ChEMBLCID, databasePassword, ChEMBLSchema):
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
        and a2t.confidence_score = '4'
        and a2t.tid = td.tid
        and td.target_type = 'PROTEIN'
        and td.db_source = 'SWISS-PROT'
    """
    # Connect to the ChEMBL schema in the database.
    conn, cursor = mysql.openConnection(databasePassword, ChEMBLSchema)
    print '\tExtracting target information.'
    cursor.execute(queryTarget2Compound)
    resultsTarget2Compound = cursor.fetchall()
    resultsTarget2Compound = list(set(resultsTarget2Compound))
    target2CompoundsDict = {}
    chemblID2Molregno = {}
    ChEMBLIDs = set([])
##    molecules = set([])
    writeOut = open(ChEMBLTargets, 'w')
    for i in resultsTarget2Compound:
##        molecules.add(i[1])
##        if target2CompoundsDict.has_key(i[0]):
##            target2CompoundsDict[i[0]].append(str(i[1]))
##        else:
##            target2CompoundsDict[i[0]] = []
##            target2CompoundsDict[i[0]].append(str(i[1]))
        chemblID2Molregno[str(i[2])] = str(i[1])
        ChEMBLIDs.add(str(i[2]))
        writeOut.write(i[0] + '\t' + str(i[1]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i[5]) + '\t' + str(i[6]) + '\t' + str(i[7]) + '\n')
##    for i in target2CompoundsDict.keys():
##        target2CompoundsDict[i] = ';'.join(set(target2CompoundsDict[i]))
##    targets = [(i + '\t' + target2CompoundsDict[i]) for i in target2CompoundsDict.keys()]
##    print '\tNumber of approved targets in ChEMBL: ', len(targets)
    writeOut.close()
##    utilities.list2file.main(resultsTarget2Compound, ChEMBLTargets)
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
    
##    return molecules

##def extract_drugs(molecules, ChEMBLDrugs, databasePassword, ChEMBLSchema):
##    # Connect to the ChEMBL schema in the database.
##    conn, cursor = mysql.openConnection(databasePassword, ChEMBLSchema)
##    print '\tExtracting drug information.'
##    writeDrugs = open(ChEMBLDrugs, 'w')
##    activities = set([])
##    print len(activities)
##    for m in molecules:
##        cursor = mysql.tableSELECT(cursor, 'pref_name', 'molecule_dictionary', 'molregno = "' + str(m) + '"')
##        prefNames = cursor.fetchall()[0][0]
##        cursor = mysql.tableSELECT(cursor, 'standard_inchi', 'compound_structures', 'molregno = "' + str(m) + '"')
##        standardInChis = cursor.fetchall()[0][0]
##        cursor = mysql.tableSELECT(cursor, 'activity_id', 'activities', 'molregno = "' + str(m) + '"')
##        activityIDs = cursor.fetchall()
##        activityIDs = [str(i[0]) for i in activityIDs]
##        activities = activities.union(activityIDs)
##        activityIDs = ';'.join(activityIDs)
##        writeDrugs.write(str(m) + '\t' + prefNames + '\t' + standardInChis + '\t' + activityIDs + '\n')
##    writeDrugs.close()
##    mysql.closeConnection(conn, cursor)
##    
##    return activities

##def extract_activities(activities, ChEMBLActivities, databasePassword, ChEMBLSchema):
##    # Connect to the ChEMBL schema in the database.
##    conn, cursor = mysql.openConnection(databasePassword, ChEMBLSchema)
##    print '\tExtracting activity information.'
##    writeActivities = open(ChEMBLActivities, 'w')
##    for i in activities:
##        cursor = mysql.tableSELECT(cursor, '*', 'activities', 'activity_id="' + i + '"')
##        resultActivities = cursor.fetchall()
##        for j in resultActivities:
##            type = str(j[11])
##            typeLower = type.lower()
##            if 'ki' in typeLower or 'kd' in typeLower or 'ic50' in typeLower or 'ec50' in typeLower:
##                # Only interested in these four affinity measures.
##                relation = str(j[5])
##                value = str(j[8])
##                units = str(j[9])
##                writeActivities.write(i + '\t' + '\t'.join([relation, value, units, type]) + '\n')
##    mysql.closeConnection(conn, cursor)