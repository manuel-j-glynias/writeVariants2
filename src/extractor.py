import os
import datetime

import mysql.connector
from src.graphql_utils import send_query
from src.sql_utils import get_local_db_connection, maybe_create_and_select_database, drop_table_if_exists, create_table, \
    get_schema, write_load_files_using_func, load_table, get_load_dir

import re


def remove_control_characters(str):
    return re.sub(r'[\x00-\x1f\x7f-\x9f]', '', str)

def get_current_author_data(server:str)->list:
    query = '{Author { id, surname, first_initial } }'

    response = send_query(query, server)
    return response['data']['Author']

def get_current_journal_data(server:str)->list:
    query = '{ Journal { id,name } }'

    response = send_query(query, server)
    return response['data']['Journal']

def get_current_literature_reference_data(server:str)->list:
    query = '{ LiteratureReference { id, PMID,DOI,title,journal{id},volume,first_page,last_page,publication_Year,shortReference,abstract,authors {id}} }'

    response = send_query(query, server)
    return response['data']['LiteratureReference']

def get_current_internet_reference_data(server:str)->list:
    query = '{ InternetReference { id, accessed_date,web_address,shortReference} }'

    response = send_query(query, server)
    return response['data']['InternetReference']


def get_current_user_data(server:str)->list:
    query = '{ User {id, name, password, isAdmin } }'

    response = send_query(query, server)
    return response['data']['User']

def get_current_editable_statement_data(server:str)->list:
    query = '{ EditableStatement { id, field, statement, edit_date,  deleted, editor { id }, references { id, __typename } } }'
    response = send_query(query, server)
    return response['data']['EditableStatement']

def get_current_jax_gene_data(server:str)->list:
    query = '{ JaxGene { id, name, canonicalTranscript, chromosome, entrezId, jaxId, synonyms, description { id } } }'
    response = send_query(query, server)
    return response['data']['JaxGene']

def get_current_myGene_info_gene_data(server:str)->list:
    query = '{ MyGeneInfo_Gene {id, name, chromosome, strand, start, end, entrezId, description { id }}}'
    response = send_query(query, server)
    return response['data']['MyGeneInfo_Gene']


def get_current_uniprot_entry_data(server:str)->list:
    query = '{ Uniprot_Entry {id, name,  accessionNumber, uniprot_id, function { id }, gene { id} }}'
    response = send_query(query, server)
    return response['data']['Uniprot_Entry']


def get_current_omnigene_data(server:str)->list:
    query = '{ OmniGene {id, name, panelName, geneDescription {id}, oncogenicCategory {id }, synonymsString {id }, myGeneInfoGene {id}, jaxGene {id}, uniprot_entry { id }} }'
    response = send_query(query, server)
    return response['data']['OmniGene']


def convert_authors(author_dict)->list:
    row = [author_dict['surname'], author_dict['first_initial'], author_dict['id']]
    return [row]

def convert_journals(j_dict)->list:
    row = [j_dict['name'],j_dict['id']]
    return [row]

def convert_references(reference_dict)->list:
    row = [reference_dict['PMID'], reference_dict['DOI'], reference_dict['title'], reference_dict['journal']['id'], reference_dict['volume'], reference_dict['first_page'],
           reference_dict['last_page'], reference_dict['publication_Year'], reference_dict['shortReference'], reference_dict['abstract'], reference_dict['id']]
    return [row]

def convert_author_refs(reference_dict)->list:
    rows = []
    for author in reference_dict['authors']:
        author_id = author['id']
        row = [None, author_id, reference_dict['id']]
        rows.append(row)
    return rows


def convert_internet_references(internet_reference_dict)->list:
    row = [internet_reference_dict['accessed_date'], internet_reference_dict['web_address'], internet_reference_dict['shortReference'], internet_reference_dict['id']]
    return [row]

def convert_users(users_dict)->list:
    isAdmin = 'TRUE'
    if users_dict['isAdmin'] == False:
        isAdmin = 'FALSE'
    row = [users_dict['name'], users_dict['password'], isAdmin, users_dict['id']]
    return [row]

def convert_editable_statements(es_dict)->list:
    row = [es_dict['field'],remove_control_characters(es_dict['statement']),es_dict['edit_date'],es_dict['editor']['id'],es_dict['id']]
    return [row]

def convert_es_refs(es_dict)->list:
    rows = []
    for ref in es_dict['references']:
        if ref['__typename'] == 'LiteratureReference':
            row = [None, es_dict['id'], ref['id']]
            rows.append(row)
    return rows

def convert_es_internet_refs(es_dict)->list:
    rows = []
    for ref in es_dict['references']:
        if ref['__typename'] == 'InternetReference':
            row = [None, es_dict['id'], ref['id']]
            rows.append(row)
    return rows

def convert_jax_genes(jax_gene_dict)->list:
    row = [jax_gene_dict['name'], jax_gene_dict['canonicalTranscript'][0], jax_gene_dict['chromosome'], jax_gene_dict['entrezId'], jax_gene_dict['jaxId'], jax_gene_dict['description']['id'], jax_gene_dict['id']]
    return [row]

def convert_jax_synonyms(jax_gene_dict)->list:
    rows = []
    for syn in jax_gene_dict['synonyms']:
        row = [None, syn, jax_gene_dict['id']]
        rows.append(row)
    return rows

def convert_myGeneInfo_genes(myGene_dict)->list:
    row = [myGene_dict['name'], myGene_dict['chromosome'], myGene_dict['strand'], myGene_dict['start'], myGene_dict['end'], myGene_dict['entrezId'], myGene_dict['description']['id'], myGene_dict['id']]
    return [row]

def convert_uniprot(uniprot_doct)->list:
    row = [uniprot_doct['name'], uniprot_doct['accessionNumber'], uniprot_doct['uniprot_id'], uniprot_doct['function']['id'], uniprot_doct['gene']['id'], uniprot_doct['id']]
    return [row]

def get_id_helper(obj,key):
    val = ''
    if key in obj:
        obj2 = obj[key]
        if obj2 != None and 'id' in obj2:
            val = obj2['id']
    return val


def convert_omnigene(omnigene_dict)->list:
    row = [omnigene_dict['name'], omnigene_dict['panelName'], get_id_helper(omnigene_dict, 'geneDescription'), get_id_helper(omnigene_dict, 'oncogenicCategory'), get_id_helper(omnigene_dict, 'synonymsString'),
           get_id_helper(omnigene_dict, 'myGeneInfoGene'), get_id_helper(omnigene_dict, 'jaxGene'), get_id_helper(omnigene_dict, 'uniprot_entry'), omnigene_dict['id']]
    return [row]


def delete_load_files(load_dir):
    for entry in os.scandir(load_dir):
        if entry.is_file:
            os.remove(entry)


def create_data_and_conversion_dictionaries(server_read):
    dict_of_data = {}
    dict_of_conversion_functions = {}

    dict_of_data['Author'] = get_current_author_data(server_read)
    dict_of_conversion_functions['Author'] = convert_authors

    dict_of_data['Journal'] = get_current_journal_data(server_read)
    dict_of_conversion_functions['Journal'] = convert_journals

    dict_of_data['LiteratureReference'] = get_current_literature_reference_data(server_read)
    dict_of_conversion_functions['LiteratureReference'] = convert_references
    # reuse refs data for join table
    dict_of_data['LiteratureReference_Author'] = dict_of_data['LiteratureReference']
    dict_of_conversion_functions['LiteratureReference_Author'] = convert_author_refs

    dict_of_data['InternetReference'] = get_current_internet_reference_data(server_read)
    dict_of_conversion_functions['InternetReference'] = convert_internet_references

    dict_of_data['User'] = get_current_user_data(server_read)
    dict_of_conversion_functions['User'] = convert_users

    dict_of_data['EditableStatement'] = get_current_editable_statement_data(server_read)
    dict_of_conversion_functions['EditableStatement'] = convert_editable_statements

    dict_of_data['EditableStatement_LiteratureReference'] = dict_of_data['EditableStatement']
    dict_of_conversion_functions['EditableStatement_LiteratureReference'] = convert_es_refs

    dict_of_data['EditableStatement_InternetReference'] = dict_of_data['EditableStatement']
    dict_of_conversion_functions['EditableStatement_InternetReference'] = convert_es_internet_refs

    dict_of_data['JaxGene'] = get_current_jax_gene_data(server_read)
    dict_of_conversion_functions['JaxGene'] = convert_jax_genes

    dict_of_data['JaxGene_Synonym'] = dict_of_data['JaxGene']
    dict_of_conversion_functions['JaxGene_Synonym'] = convert_jax_synonyms

    dict_of_data['MyGeneInfoGene'] = get_current_myGene_info_gene_data(server_read)
    dict_of_conversion_functions['MyGeneInfoGene'] = convert_myGeneInfo_genes

    dict_of_data['UniprotEntry'] = get_current_uniprot_entry_data(server_read)
    dict_of_conversion_functions['UniprotEntry'] = convert_uniprot

    dict_of_data['OmniGene'] = get_current_omnigene_data(server_read)
    dict_of_conversion_functions['OmniGene'] = convert_omnigene

    return dict_of_conversion_functions, dict_of_data

def main():
    print(datetime.datetime.now().strftime("%H:%M:%S"))
    my_db = None
    my_cursor = None
    server_read: str = '165.227.89.140'
    descriptions_csv_path = '../config/table_descriptions_02_01.csv'
    db_dict = get_schema(descriptions_csv_path)
    load_dir = get_load_dir() + 'extracted/'

    dict_of_conversion_functions, dict_of_data = create_data_and_conversion_dictionaries(server_read)
    delete_load_files(load_dir)
    write_load_files_using_func(db_dict, dict_of_data, dict_of_conversion_functions,load_dir)
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
    print(datetime.datetime.now().strftime("%H:%M:%S"))





if __name__ == "__main__":
    main()