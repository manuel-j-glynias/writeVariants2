import csv
import datetime
import mysql.connector
from src.graphql_utils import erase_neo4j, send_schema_request, send_mutation
from src.sql_utils import get_schema, get_local_db_connection


def get_mutations(server):
    query = '''{
  __schema {
    mutationType {
      name
      fields {
        name
        args {
          name
          type {
            kind
            ofType {
              kind
              name
              description
              ofType {
                name
                kind
                ofType{
                  name
                }
              }
            }
          }
        }
      }
    }
  }
}
'''
    response = send_schema_request(query, server)
    fields = response['data']['__schema']['mutationType']['fields']
    mutations = {}
    for m in fields:
        name = m['name']
        mutations[name] = m['args']
    return mutations


def get_col_head(col: str):
    l = col.split('_')
    head = l[0].lower()
    return head


def add_list_to_graphql_objects(my_cursor, mutations, method, db_dict, db_name, table_name, server):
    args = mutations[method]
    method_start = method
    if method_start.startswith('add'):
        method_start = method_start[3:]

    args_dict = {}
    for arg in args:
        add_arg_to_dict(arg, args_dict)
    list_arg = get_list_arg(args_dict)
    obj_arg = get_obj_arg(args_dict)

    cols = db_dict[db_name][table_name]['col_order']
    query = ''
    order_by = None
    object_id_index = 0
    added_id_index = 0
    i = 0
    for col in cols:
        field_info = db_dict[db_name][table_name][col]
        # item 3 is auto-incrmenet
        if field_info[3] == 'N':
            if query != '':
                query += ','
            query += col
            # col_head = get_col_head(col)
            if col.endswith('_graph_id'):
                object_name = col[:-9]
                if method_start.startswith(object_name):
                    object_id_index = i
                    order_by = col
                else:
                    added_id_index = i
            i += 1
    query = 'SELECT ' + query + ' FROM ' + db_name + '.' + table_name
    if order_by:
        query += ' ORDER BY ' + order_by
    my_cursor.execute(query)
    rows = my_cursor.fetchall()
    object_id = rows[0][object_id_index]
    counter = 0
    m = ''
    list_str = '['
    print(object_id)
    rows.append([0, 0])
    for row in rows:
        if row[object_id_index] != object_id:
            list_str += ']'
            s = table_name + '_' + str(counter) + ': ' + method + '(' + obj_arg + ': \\"' + object_id + '\\", ' + list_arg + ': ' + list_str + '), '
            m += s
            counter += 1
            if (counter % 100 == 0):
                print(m)
                send_mutation(m, server)
                m = ''
            list_str = '['
            object_id = row[object_id_index]
        list_str += f'\\"{row[added_id_index]}\\",'


def get_list_arg(args_dict):
    list_arg = None
    for arg in args_dict.values():
        if arg['is_list']:
            list_arg = arg['name']
            break
    return list_arg



def get_obj_arg(args_dict):
    obj_arg = None
    for arg in args_dict.values():
        if not arg['is_list']:
            obj_arg = arg['name']
            break
    return obj_arg


def add_to_graphql_objects(my_cursor, mutations, method, db_dict, db_name, table_name, server):
    args = mutations[method]
    args_dict = {}
    for arg in args:
        add_arg_to_dict(arg, args_dict)
    list_arg = get_list_arg(args_dict)
    obj_arg = get_obj_arg(args_dict)

    cols = db_dict[db_name][table_name]['col_order']
    query = ''
    for col in cols:
        field_info = db_dict[db_name][table_name][col]
        # item 3 is auto-incrmenet
        if field_info[3] == 'N':
            if col.lower().startswith(list_arg.lower()):
                query += col
    query += ', graph_id'
    query = 'SELECT ' + query + ' FROM ' + db_name + '.' + table_name
    counter = 0
    m = ''

    my_cursor.execute(query)
    rows = my_cursor.fetchall()
    for row in rows:
        adder_id = row[0]
        object_id = row[1]
        s = table_name + '_' + str(counter) + ': ' + method + '(' + obj_arg + ': \\"' + object_id + '\\", ' + list_arg + ': [\\"' + adder_id + '\\"]), '
        m += s
        counter += 1
        if (counter % 100 == 0):
            print(m)
            send_mutation(m, server)
            m = ''
    if m != '':
        print(m)
        send_mutation(m, server)


def create_graphql_objects(my_cursor, mutations, method, db_dict, db_name, table_name, server):
    cols = db_dict[db_name][table_name]['col_order']
    query = ''
    for col in cols:
        if query != '':
            query += ','
        query += col
    query = 'SELECT ' + query + ' FROM ' + db_name + '.' + table_name
    my_cursor.execute(query)
    rows = my_cursor.fetchall()
    args = mutations[method]
    args_dict = {}
    for arg in args:
        add_arg_to_dict(arg, args_dict)
    counter = 0
    m = ''
    for row in rows:
        s = table_name + '_' + str(counter) + ': ' + method + '('
        counter += 1
        for i in range(len(cols)):
            col = cols[i]
            val = row[i]
            if col == 'graph_id':
                col = 'id'
            if col in args_dict:
                arg = args_dict[col]
                if arg['is_boolean']:
                    val = str(val).lower()
                if arg['is_enum']:
                    val = str(val).capitalize()
                if arg['is_string']:
                    s += col + ': \\"' + str(val) + '\\", '
                else:
                    s += col + ': ' + str(val) + ', '
        s += '),'
        m += s
        if (counter % 100 == 0):
            print(m)
            send_mutation(m, server)
            m = ''
    if m != '':
        print(m)
        send_mutation(m, server)

def create_graphql_objects_with_list(my_cursor, mutations, method, db_dict, db_name, table_name, syn_dict,server):
    cols = db_dict[db_name][table_name]['col_order']
    query = ''
    for col in cols:
        if query != '':
            query += ','
        query += col
    query = 'SELECT ' + query + ' FROM ' + db_name + '.' + table_name
    my_cursor.execute(query)
    rows = my_cursor.fetchall()
    args = mutations[method]
    args_dict = {}
    for arg in args:
        add_arg_to_dict(arg, args_dict)
    counter = 0
    m = ''
    for row in rows:
        s = table_name + '_' + str(counter) + ': ' + method + '('
        counter += 1
        for i in range(len(cols)):
            col = cols[i]
            val = row[i]
            val_list = []
            if col == 'graph_id':
                col = 'id'
                graph_id = val
                if graph_id in syn_dict:
                    val_list = syn_dict[graph_id]
            if col in args_dict:
                arg = args_dict[col]
                if arg['is_boolean']:
                    val = str(val).lower()
                if arg['is_string']:
                    s += col + ': \\"' + str(val) + '\\", '
                else:
                    s += col + ': ' + str(val) + ', '
        s += 'list: ['
        for v in val_list:
            s += '\\"' + v + '\\",'
        s += ']'
        s += '),'
        m += s
        if (counter % 100 == 0):
            print(m)
            send_mutation(m, server)
            m = ''
    if m != '':
        print(m)
        send_mutation(m, server)


def add_arg_to_dict(arg, args_dict):
    name = arg['name']
    is_list = False
    if 'ofType' in arg['type'] and type(arg['type']['ofType']) is dict and 'kind' in arg['type']['ofType']:
        is_list = arg['type']['ofType']['kind'] == 'LIST'
    is_enum = 'ofType' in arg['type'] and type(arg['type']['ofType']) is dict and 'kind' in arg['type']['ofType'] and arg['type']['ofType']['kind'] == 'ENUM'
    is_int = 'ofType' in arg['type'] and type(arg['type']['ofType']) is dict and 'name' in arg['type']['ofType'] and arg['type']['ofType']['name'] == 'Int'
    is_boolean = 'ofType' in arg['type'] and type(arg['type']['ofType']) is dict and 'name' in arg['type']['ofType'] and arg['type']['ofType']['name'] == 'Boolean'
    is_string = not (is_enum or is_int or is_boolean)
    args_dict[name] = {'is_list': is_list, 'is_string': is_string, 'is_enum': is_enum, 'is_int': is_int, 'is_boolean': is_boolean, 'name':name}


def get_syn_dict(path):
    syn_dict = {}
    firstline = True
    with open(path) as csvfile:
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

    return syn_dict


def write_graphql():
    print(datetime.datetime.now().strftime("%H:%M:%S"))

    schema_graphql = '../config/schema.graphql'
    server_write: str = 'localhost'
    # server_write:str ='161.35.115.213'

    my_db = None
    my_cursor = None
    erase_neo4j(schema_graphql, server_write)
    mutations = get_mutations(server_write)
    db_dict = get_schema('../config/variant_table_descriptions.csv')

    try:
        my_db = get_local_db_connection()
        my_cursor = my_db.cursor(buffered=True)

        create_graphql_objects(my_cursor, mutations, 'createUser', db_dict, 'OmniSeqKnowledgebase2', 'User', server_write)
        create_graphql_objects(my_cursor, mutations, 'createAuthor', db_dict, 'OmniSeqKnowledgebase2', 'Author', server_write)
        create_graphql_objects(my_cursor, mutations, 'createJournal', db_dict, 'OmniSeqKnowledgebase2', 'Journal', server_write)
        create_graphql_objects(my_cursor, mutations, 'createLiteratureReference', db_dict, 'OmniSeqKnowledgebase2', 'LiteratureReference', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addLiteratureReferenceJournal', db_dict, 'OmniSeqKnowledgebase2', 'LiteratureReference', server_write)
        add_list_to_graphql_objects(my_cursor, mutations, 'addLiteratureReferenceAuthors', db_dict, 'OmniSeqKnowledgebase2', 'LiteratureReference_Author', server_write)
        create_graphql_objects(my_cursor, mutations, 'createInternetReference', db_dict, 'OmniSeqKnowledgebase2', 'InternetReference', server_write)
        create_graphql_objects(my_cursor, mutations, 'createEditableStatement', db_dict, 'OmniSeqKnowledgebase2', 'EditableStatement', server_write)
        add_list_to_graphql_objects(my_cursor, mutations, 'addEditableStatementReferences', db_dict, 'OmniSeqKnowledgebase2', 'EditableStatement_LiteratureReference', server_write)
        add_list_to_graphql_objects(my_cursor, mutations, 'addEditableStatementReferences', db_dict, 'OmniSeqKnowledgebase2', 'EditableStatement_InternetReference', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addEditableStatementEditor', db_dict, 'OmniSeqKnowledgebase2', 'EditableStatement', server_write)
        syn_dict = get_syn_dict('../load_files/Synonym.csv')
        create_graphql_objects_with_list(my_cursor, mutations, 'createEditableSynonymList', db_dict, 'OmniSeqKnowledgebase2', 'EditableSynonymList', syn_dict,server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addEditableSynonymListEditor', db_dict, 'OmniSeqKnowledgebase2', 'EditableSynonymList', server_write)
        create_graphql_objects(my_cursor, mutations, 'createJaxGene', db_dict, 'OmniSeqKnowledgebase2', 'JaxGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addJaxGeneCanonicalTranscript', db_dict, 'OmniSeqKnowledgebase2', 'JaxGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addJaxGeneDescription', db_dict, 'OmniSeqKnowledgebase2', 'JaxGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addJaxGeneSynonyms', db_dict, 'OmniSeqKnowledgebase2', 'JaxGene', server_write)

        create_graphql_objects(my_cursor, mutations, 'createMyGeneInfoGene', db_dict, 'OmniSeqKnowledgebase2', 'MyGeneInfoGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addMyGeneInfoGeneDescription', db_dict, 'OmniSeqKnowledgebase2', 'MyGeneInfoGene', server_write)
        create_graphql_objects(my_cursor, mutations, 'createUniprotEntry', db_dict, 'OmniSeqKnowledgebase2', 'UniprotEntry', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addUniprotEntryFunction', db_dict, 'OmniSeqKnowledgebase2', 'UniprotEntry', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addUniprotEntryGene', db_dict, 'OmniSeqKnowledgebase2', 'UniprotEntry', server_write)

        create_graphql_objects(my_cursor, mutations, 'createOmniGene', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        # add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneTranscript', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneGeneDescription', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneMyGeneInfoGene', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        # add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneJaxGene', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        # add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneUniprotEntry', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        # add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneOncogenicCategory', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)
        # add_to_graphql_objects(my_cursor, mutations, 'addOmniGeneSynonyms', db_dict, 'OmniSeqKnowledgebase2', 'OmniGene', server_write)









    except mysql.connector.Error as error:
        print("Failed in MySQL: {}".format(error))
    finally:
        if (my_db.is_connected()):
            my_cursor.close()
    print(datetime.datetime.now().strftime("%H:%M:%S"))


if __name__ == "__main__":
    write_graphql()
