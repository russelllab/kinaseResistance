#!/usr/bin/env python3

'''
Script to make a database for
kinase mutations project
'''

import mysql.connector
import os

def connection():
    '''Function to connect to the database'''
    mydb = mysql.connector.connect(
    host="localhost",
    user="kinase_user",
    password=""
    )
    return mydb

def create_ptm_table(mycursor)->None:
    '''Function to create the PTM table'''
    mycursor.execute("DROP TABLE IF EXISTS ptm")
    mycursor.execute("CREATE TABLE ptm (id INT AUTO_INCREMENT PRIMARY KEY, acc VARCHAR(10), gene VARCHAR(10), wtAA VARCHAR(1), wtPos INT, ptm_type VARCHAR(5), hmm_pos INT, UNIQUE(acc, gene, wtAA, wtPos, ptm_type, hmm_pos))")
    for line in open('../data/Kinase_psites_hits_split_trimmed.tsv', 'r'):
        if line.startswith('#'): continue
        line = line.rstrip().split('\t')
        acc = line[0]
        gene = line[1]
        ptm_type = line[3].split('-')[1]
        residue = line[3].split('-')[0]
        wtPos = int(residue[1:])
        wtAA = residue[0]
        hmm_pos = line[4]
        mycursor.execute("INSERT INTO ptm (acc, gene, wtAA, wtPos, ptm_type, hmm_pos) VALUES (%s, %s, %s, %s, %s, %s)", (acc, gene, wtAA, wtPos, ptm_type, hmm_pos))

if __name__ == '__main__':
    mydb = connection()
    mycursor = mydb.cursor(buffered=True)
    # Execute SQL query to check if database exists
    mycursor.execute("SHOW DATABASES")

    # Loop through results and check if database exists
    db_exists = False
    for db in mycursor:
        if db[0] == "mydatabase":
            db_exists = True
            break
    if db_exists == False: mycursor.execute("CREATE DATABASE mydatabase")
    mycursor.execute("use mydatabase")

    # Create PTM table
    create_ptm_table(mycursor)
    mydb.commit()

    # Use mysqldump to create backup file
    backup_file = "kinaseDB.sql"
    os.system(f"mysqldump -u {mydb.user} {mydb.database} > {backup_file}")

    # Close MySQL connection
    mydb.close()
    # main()