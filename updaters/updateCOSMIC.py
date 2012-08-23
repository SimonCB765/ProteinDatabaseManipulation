import utilities.file2list
import utilities.MySQLaccess as mysql

def main(cosmicParsedGene, cosmicParsedGene2Mutation, cosmicParsedMutation, databasePassword, schemaProteins,
         tableCOSMICGene, tableCOSMICGene2Mutation, tableCOSMICMutation):

    #===========================================================================
    # Extract and format the parsed data.
    #===========================================================================
    cosmicGene = utilities.file2list.main(cosmicParsedGene)
    cosmicGene = [eval(i) for i in cosmicGene]
    cosmicGene2Mutation = utilities.file2list.main(cosmicParsedGene2Mutation)
    cosmicGene2Mutation = [eval(i) for i in cosmicGene2Mutation]
    cosmicMutation = utilities.file2list.main(cosmicParsedMutation)
    cosmicMutation = [eval(i) for i in cosmicMutation]

    #===========================================================================
    # Update the database tables.
    #===========================================================================
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database=schemaProteins)
    values = '(' + ('%s,' * len(cosmicGene[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableCOSMICGene, values, cosmicGene)

    values = '(' + ('%s,' * len(cosmicMutation[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableCOSMICMutation, values, cosmicMutation)

    values = '(' + ('%s,' * len(cosmicGene2Mutation[0]))
    values = values[:-1] + ')'
    cursor = mysql.tableINSERT(cursor, tableCOSMICGene2Mutation, values, cosmicGene2Mutation)
    mysql.closeConnection(conn, cursor)