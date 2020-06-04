import csv
import shutil
from typing import Callable, List, Dict
from pathlib import Path
import mysql.connector
import sys
import os


def get_sql_credentials():
    in_file = open('../config/sql.txt', 'r')
    for line in in_file:
        item = line.split('=')
        if (item[0] == 'user'):
            sql_user = item[-1].strip()
        elif (item[0] == 'passwd'):
            sql_passwd = item[-1].strip()
    return(sql_user,sql_passwd)

def get_cloud_db_connection():
    sql_user, sql_passwd = get_sql_credentials()
    my_db = mysql.connector.connect(
        host="159.203.79.49",
        user=sql_user,
        passwd=sql_passwd
    )
    return my_db


def get_local_db_connection():
    sql_user, sql_passwd = get_sql_credentials()
    my_db = mysql.connector.connect(
        host="localhost",
        user=sql_user,
        passwd=sql_passwd,
        allow_local_infile = True
    )
    return my_db


def drop_database(my_cursor,db_name):
    if does_db_exist(my_cursor,db_name):
        my_cursor.execute('DROP DATABASE ' + db_name)


def maybe_create_and_select_database(my_cursor, db_name):
    if not does_db_exist(my_cursor, db_name):
        my_cursor.execute('CREATE DATABASE ' + db_name + ' DEFAULT CHARACTER SET utf8mb4')
        print(f"{db_name} Database created successfully")
    my_cursor.execute('USE ' + db_name)


def does_db_exist(my_cursor,db_name):
    exists = False
    my_cursor.execute("SHOW DATABASES")
    for x in my_cursor:
        if x[0].lower() == db_name.lower():
            exists = True
            break
    return exists


def does_table_exist(my_cursor,table_name):
    exists = False
    my_cursor.execute("SHOW TABLES")
    for x in my_cursor:
        if x[0] == table_name:
            exists = True
            break
    return exists

def drop_table_if_exists(my_cursor,table_name):
    sql = "DROP TABLE IF EXISTS " + table_name
    my_cursor.execute(sql)

def create_table(my_cursor,table_name,db_name,db_dict):
    sql = 'CREATE TABLE ' + table_name + ' ('
    primary_keys = []
    indexes = []
    for col in db_dict[db_name][table_name]['col_order']:
        # determine data type for column
        col_type = db_dict[db_name][table_name][col][0]

        # check if a primary key or an index
        col_key = db_dict[db_name][table_name][col][1]
        if (col_key == 'Primary'):
            primary_keys.append(col)
        if (col_key == 'Index'):
            indexes.append(col)

        # check if NULL values allowed
        if (db_dict[db_name][table_name][col][2] == 'N'):
            null_stmt = ' NOT NULL '
        else:
            null_stmt = ''

        # check if field should AUTO_INCREMENT
        if (db_dict[db_name][table_name][col][3] == 'Y'):
            auto_incr_stmt = ' AUTO_INCREMENT '
        else:
            auto_incr_stmt = ''

        if (col != db_dict[db_name][table_name]['col_order'][-1]):
            sql += col + ' ' + col_type + null_stmt + auto_incr_stmt + ', '
        else:
            sql += col + ' ' + col_type + null_stmt + auto_incr_stmt

    if (len(primary_keys) > 0):
        sql += ', PRIMARY KEY (' + ', '.join(primary_keys) + ')'

    if (len(indexes) > 0):
        for cur_index_col in indexes:
            sql += ', INDEX ' + cur_index_col + '_idx (' + cur_index_col + ')'
    sql += ')'

    my_cursor.execute(sql)

def load_table(my_cursor, table_name, col_list, load_dir):
    load_path = load_dir + table_name + '.csv'
    num_cols = len(col_list)

    sql0 = 'TRUNCATE TABLE ' + table_name
    sql1 = make_sql('utf8mb4', load_path, table_name, col_list)
    print('Loading', table_name)

    my_cursor.execute(sql0)
    try:
        my_cursor.execute(sql1)
    except Exception as e:
        # if ('Incorrect string value' in e.args[1]):
        if ((e.args[0] == 1300) | (e.args[0] == 1366)):
            print('WARNING: Load failed using utfmb4 charset, try with cp1251...')
            sql1_alt = make_sql('cp1251', load_path, table_name, col_list)
            try:
                my_cursor.execute(sql1_alt)
                print('SUCCESS!!!')
            except Exception as e:
                print(e)
                print(sql1_alt)
                print('ERROR: Load failed using cp1251 charset...')
                sys.exit()
        else:
            print('ERROR:', e)
            print(sql1)
            sys.exit()


# load_sql = "LOAD DATA LOCAL INFILE '/tmp/city.csv' INTO TABLE usermanaged.city\
#  FIELDS TERMINATED BY ',' ENCLOSED BY '"' IGNORE 1 LINES;"
def make_sql(char_set, load_path, table_name, col_list):
    sql = 'LOAD DATA LOCAL INFILE \'' + load_path + '\' INTO TABLE ' + table_name
    sql += ' CHARACTER SET ' + char_set + ' FIELDS TERMINATED BY \',\' OPTIONALLY ENCLOSED BY \'\"\'  IGNORE 1 LINES ('
    # This will set entry to NULL if encounter a blank field in load file
    for col in col_list:
        if (col != col_list[-1]):
            sql += '@' + col + ', '
        else:
            sql += '@' + col + ') SET '

    for col in col_list:
        if (col != col_list[-1]):
            sql += col + ' = nullif(@' + col + ',\'\'), '
        else:
            sql += col + ' = nullif(@' + col + ',\'\') '

    return sql

def get_schema(descriptions_csv_path):
    db_dict = {}  # {db_name:{table_name:{col:[type,key,allow_null,ref_col_list],'col_order':[cols in order]}}}

    init_file = open(descriptions_csv_path, 'r')
    reader = csv.reader(init_file, quotechar='\"')
    for line in reader:
        #print(line)
        db_name = line[0]
        if (db_name == 'Database'):
            continue
        if (db_name not in db_dict.keys()):
            db_dict[db_name] = {}
        table_name = line[1]
        col = line[2]
        col_type = line[3]
        col_key = line[4]
        allow_null = line[5]
        auto_incr = line[6]
        ref_col_list = line[7].split('|') # we will ignore this for now during development
        try:
            ref_col_list.remove('')
        except:
            pass

        try:
            db_dict[db_name][table_name][col] = [col_type, col_key, allow_null, auto_incr, ref_col_list]
            db_dict[db_name][table_name]['col_order'].append(col)
        except:
            db_dict[db_name][table_name] = {col: [col_type, col_key, allow_null, auto_incr, ref_col_list]}
            db_dict[db_name][table_name] = {col: [col_type, col_key, allow_null, auto_incr, ref_col_list], 'col_order': [col]}
    init_file.close()
    #pprint(db_dict)
    return db_dict

def write_load_files_using_func(db_dict,dict_of_data,dict_of_conversion_functions, load_dir):

    for db_name in sorted(db_dict.keys()):
        for table_name in sorted(db_dict[db_name].keys()):
            out_file = open(load_dir + table_name + '.csv', 'w', encoding='utf-8')
            header = db_dict[db_name][table_name]['col_order']
            writer = csv.writer(out_file, lineterminator='\n')
            writer.writerow(header)
            if table_name in dict_of_conversion_functions and table_name in dict_of_data:
                func = dict_of_conversion_functions[table_name]
                data = dict_of_data[table_name]
                for row in data:
                    data_rows = func(row)
                    for r in data_rows:
                        writer.writerow(r)
            out_file.close()


def get_load_dir():
    getcwd = os.getcwd()
    if (('\\src' in getcwd) | ('/src' in getcwd)):  # this handles Windows or *nix paths
        home_dir = getcwd.split('src')[0].replace('\\', '/')
        load_dir = home_dir + 'load_files/'
    else:
        print('Need to run script from src/ directory.')
        sys.exit()
    return load_dir