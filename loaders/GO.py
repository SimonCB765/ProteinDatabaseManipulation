'''
Created on 7 Oct 2011

@author: Simon Bull
'''

import subprocess

import utilities.MySQLaccess as mysql

def main(databasePassword, schemaGO, completeGODatabase, MySQLBin):
    # Connect to the MySQL database
    conn, cursor = mysql.openConnection(inputPass = databasePassword, database = 'do not connect')
    # Setup the schemas and tables needed to record the GO information. Delete schemas and tables that
    # already exist.
    mysql.schemaDrop(cursor, 'IF EXISTS ' + schemaGO)
    mysql.schemaCreate(cursor, 'IF NOT EXISTS ' + schemaGO)
    mysql.closeConnection(conn, cursor)
    # Load the downloaded GO database into the newly created GO schema.
    subprocess.call('mysql.exe -u root -p' + databasePassword + ' ' + schemaGO + ' < ' + completeGODatabase, shell=True, cwd=MySQLBin)