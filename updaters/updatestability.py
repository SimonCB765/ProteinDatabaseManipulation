'''
Created on 5 Dec 2011

@author: Simon Bull
'''

from Bio.SeqUtils.ProtParam import ProteinAnalysis

import utilities.MySQLaccess as mysql

def main(databasePassword, schemaProteins, tableProteinInfo, tableStability):

    # Define N-terminus half life values (explanation http://en.wikipedia.org/wiki/N-end_rule and the ProtParam tool).
    halfLife = {'A' : 4.4, 'C' : 1.2, 'D' : 1.1, 'E' : 1.0, 'F' : 1.1, 'G' : 30.0, 'H' : 3.5, 'I' : 20.0, 'K' : 1.3,
                'L' : 5.5, 'M' : 30.0, 'N' : 1.4, 'P' : 20.0, 'Q' : 0.8, 'R' : 1.0, 'S' : 1.9, 'T' : 7.2,
                'V' : 100.0, 'W' : 2.8, 'Y' : 2.8}

    # Extract all the sequences stored in the database.
    conn, cursor = mysql.openConnection(databasePassword, schemaProteins)
    cursor = mysql.tableSELECT(cursor, 'UPAccession, Sequence', tableProteinInfo)
    results = cursor.fetchall()

    # Calculate the half life and instability index for each protein.
    stabilityTuples = []
    for i in results:
        sequence = i[1]
        if halfLife.has_key(sequence[0]):
            protHalfLife = halfLife[sequence[0]]
        else:
            # This will occur when the N-terminal is not an amino acid with an associated half-life value (e.g. X, B, etc.)
            protHalfLife = -1
        analysedSeq = ProteinAnalysis(sequence)
        try:
            instabilityIndex = analysedSeq.instability_index()
        except:
            instabilityIndex = -1
            print '\tContains invalid aa code: ', i[0]
        stabilityTuples.append(tuple([i[0], protHalfLife, instabilityIndex]))

    cursor.execute('TRUNCATE TABLE ' + tableStability)
    values = '(' + ('%s,' * len(stabilityTuples[0]))
    values = values[:-1] + ')'
    mysql.tableINSERT(cursor, tableStability, values, stabilityTuples)
    mysql.closeConnection(conn, cursor)

#def instability_index(prot, sequence):
#
#    # A two dimentional dictionary for calculating the instability index.
#    # Guruprasad K., Reddy B.V.B., Pandit M.W.    Protein Engineering 4:155-161(1990).
#    # It is based on dipeptide values therefore the vale for the dipeptide DG is DIWV['D']['G'].
#    DIWV = {'A': {'A': 1.0, 'C': 44.94, 'E': 1.0, 'D': -7.49,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': -7.49,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
#                  'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
#            'C': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 20.26,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 33.60,
#                  'K': 1.0, 'M': 33.60, 'L': 20.26, 'N': 1.0,
#                  'Q': -6.54, 'P': 20.26, 'S': 1.0, 'R': 1.0,
#                  'T': 33.60, 'W': 24.68, 'V': -6.54, 'Y': 1.0},
#            'E': {'A': 1.0, 'C': 44.94, 'E': 33.60, 'D': 20.26,
#                  'G': 1.0, 'F': 1.0, 'I': 20.26, 'H': -6.54,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': 1.0,
#                  'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0},
#            'D': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
#                  'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0,
#                  'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 1.0, 'P': 1.0, 'S': 20.26, 'R': -6.54,
#                  'T': -14.03, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
#            'F': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 13.34,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
#                  'K': -14.03, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
#                  'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 33.601},
#            'I': {'A': 1.0, 'C': 1.0, 'E': 44.94, 'D': 1.0,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 13.34,
#                  'K': -7.49, 'M': 1.0, 'L': 20.26, 'N': 1.0,
#                  'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0,
#                  'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
#            'G': {'A': -7.49, 'C': 1.0, 'E': -6.54, 'D': 1.0,
#                  'G': 13.34, 'F': 1.0, 'I': -7.49, 'H': 1.0,
#                  'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': -7.49,
#                  'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0,
#                  'T': -7.49, 'W': 13.34, 'V': 1.0, 'Y': -7.49},
#            'H': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
#                  'G': -9.37, 'F': -9.37, 'I': 44.94, 'H': 1.0,
#                  'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 24.68,
#                  'Q': 1.0, 'P': -1.88, 'S': 1.0, 'R': 1.0,
#                  'T': -6.54, 'W': -1.88, 'V': 1.0, 'Y': 44.94},
#            'K': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
#                  'G': -7.49, 'F': 1.0, 'I': -7.49, 'H': 1.0,
#                  'K': 1.0, 'M': 33.60, 'L': -7.49, 'N': 1.0,
#                  'Q': 24.64, 'P': -6.54, 'S': 1.0, 'R': 33.60,
#                  'T': 1.0, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
#            'M': {'A': 13.34, 'C': 1.0, 'E': 1.0, 'D': 1.0,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 58.28,
#                  'K': 1.0, 'M': -1.88, 'L': 1.0, 'N': 1.0,
#                  'Q': -6.54, 'P': 44.94, 'S': 44.94, 'R': -6.54,
#                  'T': -1.88, 'W': 1.0, 'V': 1.0, 'Y': 24.68},
#            'L': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0,
#                  'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
#                  'K': -7.49, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 33.60, 'P': 20.26, 'S': 1.0, 'R': 20.26,
#                  'T': 1.0, 'W': 24.68, 'V': 1.0, 'Y': 1.0},
#            'N': {'A': 1.0, 'C': -1.88, 'E': 1.0, 'D': 1.0,
#                  'G': -14.03, 'F': -14.03, 'I': 44.94, 'H': 1.0,
#                  'K': 24.68, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': -6.54, 'P': -1.88, 'S': 1.0, 'R': 1.0,
#                  'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 1.0},
#            'Q': {'A': 1.0, 'C': -6.54, 'E': 20.26, 'D': 20.26,
#                  'G': 1.0, 'F': -6.54, 'I': 1.0, 'H': 1.0,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0,
#                  'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 1.0,
#                  'T': 1.0, 'W': 1.0, 'V': -6.54, 'Y': -6.54},
#            'P': {'A': 20.26, 'C': -6.54, 'E': 18.38, 'D': -6.54,
#                  'G': 1.0, 'F': 20.26, 'I': 1.0, 'H': 1.0,
#                  'K': 1.0, 'M': -6.54, 'L': 1.0, 'N': 1.0,
#                  'Q': 20.26, 'P': 20.26, 'S': 20.26, 'R': -6.54,
#                  'T': 1.0, 'W': -1.88, 'V': 20.26, 'Y': 1.0},
#            'S': {'A': 1.0, 'C': 33.60, 'E': 20.26, 'D': 1.0, 'G': 1.0, 'F': 1.0, 'I': 1.0, 'H': 1.0,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 20.26, 'P': 44.94, 'S': 20.26, 'R': 20.26,
#                  'T': 1.0, 'W': 1.0, 'V': 1.0, 'Y': 1.0},
#            'R': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 20.26,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': 13.34, 'Q': 20.26, 'P': 20.26, 'S': 44.94, 'R': 58.28,
#                  'T': 1.0, 'W': 58.28, 'V': 1.0, 'Y': -6.54},
#            'T': {'A': 1.0, 'C': 1.0, 'E': 20.26, 'D': 1.0, 'G': -7.49, 'F': 13.34, 'I': 1.0, 'H': 1.0,
#                  'K': 1.0, 'M': 1.0, 'L': 1.0, 'N': -14.03, 'Q': -6.54, 'P': 1.0, 'S': 1.0, 'R': 1.0,
#                  'T': 1.0, 'W': -14.03, 'V': 1.0, 'Y': 1.0},
#            'W': {'A': -14.03, 'C': 1.0, 'E': 1.0, 'D': 1.0, 'G': -9.37, 'F': 1.0, 'I': 1.0, 'H': 24.68,
#                  'K': 1.0, 'M': 24.68, 'L': 13.34, 'N': 13.34, 'Q': 1.0, 'P': 1.0, 'S': 1.0, 'R': 1.0,
#                  'T': -14.03, 'W': 1.0, 'V': -7.49, 'Y': 1.0},
#            'V': {'A': 1.0, 'C': 1.0, 'E': 1.0, 'D': -14.03, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 1.0,
#                  'K': -1.88, 'M': 1.0, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 20.26, 'S': 1.0, 'R': 1.0,
#                  'T': -7.49, 'W': 1.0, 'V': 1.0, 'Y': -6.54},
#            'Y': {'A': 24.68, 'C': 1.0, 'E': -6.54, 'D': 24.68, 'G': -7.49, 'F': 1.0, 'I': 1.0, 'H': 13.34,
#                  'K': 1.0, 'M': 44.94, 'L': 1.0, 'N': 1.0, 'Q': 1.0, 'P': 13.34, 'S': 1.0, 'R': -15.91,
#                  'T': -7.49, 'W': -9.37, 'V': 1.0, 'Y': 13.34},
#            }
#
#    score = 0.0
#    for i in range(len(sequence) - 1):
#        if DIWV.has_key(sequence[i]):
#            if DIWV[sequence[i]].has_key(sequence[i+1]):
#                score += DIWV[sequence[i]][sequence[i+1]]
#    return (10.0 / len(sequence)) * score