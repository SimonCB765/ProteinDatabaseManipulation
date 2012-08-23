'''
Created on 7 Oct 2011

@author: Simon Bull
'''

import subprocess

import utilities.MySQLaccess as mysql

def main(databasePassword, schemaChEMBL, completeChEMBLDatabase, MySQLBin):
    # Create the ChEMBL schema overwriting any schema already in the database.
    conn, cursor = mysql.openConnection(inputPass=databasePassword, database='do not connect')
    mysql.schemaDrop(cursor, 'IF EXISTS ' + schemaChEMBL)
    mysql.schemaCreate(cursor, 'IF NOT EXISTS ' + schemaChEMBL)
    mysql.closeConnection(conn, cursor)
    # Load the downloaded ChEMBL database into the newly created ChEMBL schema.
    subprocess.call('mysql.exe -u root -p' + databasePassword + ' ' + schemaChEMBL + ' < ' + completeChEMBLDatabase, shell=True, cwd=MySQLBin)