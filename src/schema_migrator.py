import csv
import datetime
import os
import shutil

import mysql.connector

from src.sql_utils import get_local_db_connection, maybe_create_and_select_database, drop_table_if_exists, create_table, \
    get_schema, write_load_files_using_func, load_table, get_load_dir


def copy_files(files_to_copy, copy_dir,destination):
    for file in files_to_copy:
        path = copy_dir + file + '.csv'
        shutil.copy(os.path.realpath(path), os.path.realpath(destination))


def create_load_files_dict(db_dict, load_dir):
    load_files_dict = {}
    for table_name in sorted(db_dict.keys()):
        out_file = open(load_dir + table_name + '.csv', 'w', encoding='utf-8')
        header = db_dict[table_name]['col_order']
        writer = csv.writer(out_file, lineterminator='\n')
        writer.writerow(header)
        load_files_dict[table_name] = {'writer':writer, 'file':out_file}
    return load_files_dict

def open_file_for_append(file_name,load_dir,load_files_dict):
    out_file = open(load_dir + file_name + '.csv', 'a', encoding='utf-8')
    writer = csv.writer(out_file, lineterminator='\n')
    load_files_dict[file_name] = {'writer': writer, 'file': out_file}

def close_load_files(load_files_dict):
    for writer in load_files_dict.values():
        writer['file'].close()

def get_loader_user_id(extract_dir):
    user_id = None
    firstline = True
    with open(extract_dir+'User.csv') as csvfile:
        users = csv.reader(csvfile)
        for row in users:
            if firstline:
                firstline = False
            else:
                if row[0]=='loader':
                    user_id = row[3]
                    break
    return user_id


def migrate_jax_gene(load_files_dict,extract_dir,user_id):
#      make dict of synonyms by jax_id
    syn_dict = {}
    firstline = True
    with open(extract_dir+'JaxGene_Synonym.csv') as csvfile:
        synonyms = csv.reader(csvfile)
        for row in synonyms:
            if firstline:
                firstline = False
            else:
                id = row[2]
                syn = row[1]
                if id not in syn_dict:
                    syn_dict[id] = []
                syn_dict[id].append(syn)
    firstline = True
    synonym_writer = load_files_dict['Synonym']['writer']
    editable_synonym_list_writer = load_files_dict['EditableSynonymList']['writer']
    jax_gene_writer = load_files_dict['JaxGene']['writer']
    editable_statement_writer = load_files_dict['EditableStatement']['writer']

    with open(extract_dir+'JaxGene.csv') as csvfile:
        jax_genes = csv.reader(csvfile)
        for row in jax_genes:
            if firstline:
                firstline = False
            else:
                now = datetime.datetime.now()
                edit_date: str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
                name = row[0]
                jax_id = row[6]
                synonyms = []
                if jax_id in syn_dict:
                    synonyms = syn_dict[jax_id]
                EditableSynonymList_graph_id = 'esl_' + now.strftime("%Y%m%d%H%M%S%f")
                for syn in synonyms:
                    synonym_writer.writerow([None,syn,EditableSynonymList_graph_id])
                esyn_field = 'synonyms_' + jax_id
                editable_synonym_list_writer.writerow([esyn_field,edit_date,user_id,EditableSynonymList_graph_id])

                canonicalTranscript = row[1]
                ct_field = 'canonicalTranscript_' + jax_id
                canonicalTranscript_graph_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
                editable_statement_writer.writerow([ct_field,canonicalTranscript,edit_date,user_id,canonicalTranscript_graph_id])

                jax_gene_writer.writerow([row[0],canonicalTranscript_graph_id,row[2],row[3],row[4],row[5],EditableSynonymList_graph_id,jax_id])


def migrate_omnigene(load_files_dict,extract_dir, user_id):
    firstline = True
    syn_dict = {}
    with open(extract_dir + 'EditableStatement.csv') as csvfile:
        synonyms = csv.reader(csvfile)
        for row in synonyms:
            if firstline:
                firstline = False
            else:
                if row[0].startswith('SynonymsString_omnigene'):
                    ss = row[1]
                    graph_id = row[4]
                    syn_dict[graph_id] = ss


    synonym_writer = load_files_dict['Synonym']['writer']
    editable_synonym_list_writer = load_files_dict['EditableSynonymList']['writer']
    omnigene_writer = load_files_dict['OmniGene']['writer']
    firstline = True
    with open(extract_dir+'OmniGene.csv') as csvfile:
        omnigenes = csv.reader(csvfile)
        for row in omnigenes:
            if firstline:
                firstline = False
            else:
                now = datetime.datetime.now()
                edit_date: str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
                omnigene_id = row[8]
                syn_es = row[4]
                synonyms = []
                # if syn_es in syn_dict:
                synonyms = syn_dict[syn_es].split(',')

                EditableSynonymList_graph_id = 'esl_' + now.strftime("%Y%m%d%H%M%S%f")
                for syn in synonyms:
                    synonym_writer.writerow([None,syn,EditableSynonymList_graph_id])
                esyn_field = 'synonyms_' + omnigene_id
                editable_synonym_list_writer.writerow([esyn_field,edit_date,user_id,EditableSynonymList_graph_id])
                omnigene_writer.writerow([row[0],row[1],row[2],row[3],EditableSynonymList_graph_id,row[5],row[6],row[7],omnigene_id])


def migrate():
    files_to_copy = ['Author','EditableStatement','EditableStatement_InternetReference',"EditableStatement_LiteratureReference",
                     'InternetReference','Journal','LiteratureReference','LiteratureReference_Author','MyGeneInfoGene',
                     'UniprotEntry','User']


    load_dir = get_load_dir()
    extract_dir = get_load_dir() + 'extracted/'
    copy_files(files_to_copy, extract_dir, load_dir)
    descriptions_csv_path = '../config/variant_table_descriptions.csv'
    db_dict = get_schema(descriptions_csv_path)['OmniSeqKnowledgebase2']

    for file in files_to_copy:
        db_dict.pop(file)
    print(db_dict)

    load_files_dict = create_load_files_dict(db_dict,load_dir)
    open_file_for_append('EditableStatement', load_dir, load_files_dict)
    user_id = get_loader_user_id(extract_dir)
    migrate_jax_gene(load_files_dict, extract_dir, user_id)
    migrate_omnigene(load_files_dict, extract_dir, user_id)

    close_load_files(load_files_dict)

    db_dict = get_schema(descriptions_csv_path)
    try:
        my_db = get_local_db_connection()
        my_cursor = my_db.cursor(buffered=True)

        for db_name in sorted(db_dict.keys()):
            maybe_create_and_select_database(my_cursor, db_name)
            for table_name in sorted(db_dict[db_name].keys()):
                drop_table_if_exists(my_cursor, table_name)
                create_table(my_cursor, table_name, db_name, db_dict)
                load_table(my_cursor, table_name, db_dict[db_name][table_name]['col_order'],load_dir)
        my_db.commit()
    except mysql.connector.Error as error:
        print("Failed in MySQL: {}".format(error))
    finally:
        if (my_db.is_connected()):
            my_cursor.close()


if __name__ == "__main__":
    migrate()