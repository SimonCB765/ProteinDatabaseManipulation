'''
Created on 6 Apr 2011

@author: Simon Bull
'''

import MySQLdb

def openConnection(inputPass = 'root', database = '-----'):
    """Opens a connection to a local database.

    Note: the local database must be set up to use TCP/IP with the default port.
    Note: you can avoid connecting to a database by setting it to 'do not connect'.
    inputPass @type - string
    inputPass @use  - the password used to access the MySQL server
    database @type - string
    database @use  - the name of the database to connect to

    return @type - a MySQLdb connection
    return @use  - the connection to the server, and optionally database, specififed
    return @type - a MySQLdb cursor
    return @use  - the cursor corresponding to the returned connection
    """

    if inputPass == '-----':
        inputPass = raw_input('Please input your root user MySQL password followed by the the enter key: ')
    if database == '-----':
        database = raw_input('Please input the database you wish to connect to followed by the the enter key: ')
        conn = MySQLdb.connect (user = 'root', passwd = inputPass, db = database)
    elif database == 'do not connect':
        conn = MySQLdb.connect (user = 'root', passwd = inputPass)
    else:
        conn = MySQLdb.connect (user = 'root', passwd = inputPass, db = database)
    cursor = conn.cursor()

    return conn, cursor

def closeConnection(conn, cursor):
    """Closes a database connection.

    Note: there is no checking for whether the cursor corresponds to the connection.
    conn @type - a MySQLdb connection
    conn @use  - a MySQLdb connection you want to close
    cursor @type - a MySQLdb cursor
    cursor @use  - a MySQLdb cursor you want to close
    """

    cursor.close()
    conn.commit()
    conn.close()

def tableSELECT(cursor, selectRow, fromTable, whereCondition = 'not needed'):
    """Perform a SELECT..FROM..WHERE query.

    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor which you wish to use to perform the query
    selectRow @type - string with no leading or trailing whitespace
    selectRow @use  - the portion of the query coming between SELECT and FROM
    fromTable @type - string with no leading or trailing whitespace
    fromTable @use  - the portion of the query coming between FROM and WHERE
    whereCondition @type - string with no leading or trailing whitespace
    whereCondition @use  - the optional portion of the query coming after WHERE

    return @type - a MySQLdb cursor
    return @use  - returns the cursor used to perform the query, and therefore the one with the results
    """

    if type(selectRow) != str or type(fromTable) != str or type(whereCondition) != str:
        raise TypeError("All inputs except the cursor must be strings.")

    if whereCondition == 'not needed':
        cursor.execute('SELECT ' + selectRow + ' FROM ' + fromTable)
    else:
        cursor.execute('SELECT ' + selectRow + ' FROM ' + fromTable + ' WHERE ' + whereCondition)

    return cursor

def tableINSERT(cursor, insertInto, values, tupleList):
    """Performs an INSERT INTO..VALUES query.

    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor which you wish to use to perform the query
    insertInto @type - string with no leading or trailing whitespace
    insertInto @use  - the portion of the query between INSERT INTO and VALUES
    values @type - string with no leading or trailing whitespace
    values @use  - the portion of the query after VALUES
    tupleList @type - a list of tuples
    tupleList @use  - each tuple in the list should be a valid tuple to insert into the database with the elements of the tuple in the order of the columns in the database

    return @type - a MySQLdb cursor
    return @use  - returns the cursor used to perform the query, and therefore the one with the results

    Example: tableINSERT(cursor, 'animal', '(%s,%s,%s)', [(45,'Elephant','N')])
    """

    if type(insertInto) != str or type(values) != str:
        raise TypeError("All inputs except the cursor and list of tuple to insert must be strings.")

    if type(tupleList) != list:
        raise TypeError("tupleList must be of type list")
    position = 0
    for i in tupleList:
        if type(i) != tuple:
            error = 'The tuple in tupleList at position ' + str(position) + ' is not a tuple'
            raise TypeError(error)
        position += 1

    cursor.executemany('INSERT INTO ' + insertInto + ' VALUES ' + values, tupleList)

    return cursor

def tableUPDATE(cursor, table, toSet, whereCon):
    """
    """

    if type(table) != str or type(toSet) != str or type(whereCon) != str:
        raise TypeError("All inputs except the cursor must be strings.")

    query = 'UPDATE ' + table + ' SET ' + toSet + ' WHERE ' + whereCon
    cursor.execute(query)

    return cursor

def tableALTER(cursor, table, alteration, column, location = 'NA'):
    """
    table is the table to add/delete to/from
    alteration is the alteration to make (i.e. add/delete)
    column is the column to use in the alteration (must include the type of values to enter in the column (i.e. VARCHAR(40))
    location is the potential keyword such as FIRST or AFTER X
    """

    query = 'ALTER TABLE ' + table + ' ' + alteration + ' ' + column

    if location != 'NA':
        query += ' '
        query += location
    cursor.execute(query)

    return cursor

def tableEXISTS(cursor, table):
    """Assumes that cursor is connected to the database from which you want to check if a table exists.
    CAN NOT give a table that has a '.' in front of it to this function"""

    cursor.execute('SHOW TABLES LIKE \'' + table + '\'')
    return cursor

def rowDELETE(cursor, table, whereCon):
    """
    """

    if type(table) != str or type(whereCon) != str:
        raise TypeError("All inputs except the cursor must be strings.")

    query = 'DELETE FROM ' + table + ' WHERE ' + whereCon
    cursor.execute(query)

    return cursor

def schemaCreate(cursor, query):
    """Performs a query to create a schema.

    Note: the cursor should not be associated with a MySQLdb connection that is connected to a specific database.
    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor you wish to use in order to add the schema
    query @type - string with no leading or trailing whitespace
    query @use  - the portion of the MySQL query that comes after 'CREATE SCHEMA '
    """

    cursor.execute('CREATE SCHEMA ' + query)

def schemaDrop(cursor, query):
    """Performs a query to drop a schema.

    Note: the cursor should not be associated with a MySQLdb connection that is connected to a specific database.
    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor you wish to use in order to drop the schema
    query @type - string with no leading or trailing whitespace
    query @use  - the portion of the MySQL query that comes after 'DROP SCHEMA '
    """

    cursor.execute('DROP SCHEMA ' + query)

def tableCreate(cursor, query):
    """Performs a query to create a table.

    Note: the cursor should be associated with a MySQLdb connection that is connected to no database or the database you wish to add the table to.
    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor you wish to use in order to add the table
    query @type - string with no leading or trailing whitespace
    query @use  - the portion of the MySQL query that comes after 'CREATE TABLE '
    """

    cursor.execute('CREATE TABLE ' + query)

def tableDrop(cursor, query):
    """Performs a query to drop a table.

    Note: the cursor should be associated with a MySQLdb connection that is connected to no database or the database from which you wish to drop the table.
    cursor @type - a MySQLdb cursor
    cursor @use  - the cursor you wish to use in order to drop the table
    query @type - string with no leading or trailing whitespace
    query @use  - the portion of the MySQL query that comes after 'DROP TABLE '
    """

    cursor.execute('DROP TABLE ' + query)
