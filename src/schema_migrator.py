import csv
import datetime
import os
import shutil
import json
import mysql.connector

from src.graphql_utils import replace_characters, get_reference_from_pmid_by_metapub, fix_author_id, get_authors_names, ref_name_from_authors_pmid_and_year
from src.sql_utils import get_local_db_connection, maybe_create_and_select_database, drop_table_if_exists, create_table, \
    get_schema, write_load_files_using_func, load_table, get_load_dir


def copy_files(files_to_copy, copy_dir,destination):
    for file in files_to_copy:
        path = copy_dir + file + '.csv'
        shutil.copy(os.path.realpath(path), os.path.realpath(destination))


def create_load_files_dict(db_dict, load_dir):
    load_files_dict = {}
    for table_name in sorted(db_dict.keys()):
        out_file = open(load_dir + table_name + '.csv', 'a', encoding='utf-8')
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

    jax_transcript_dict = {}
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
                jax_transcript_dict[jax_id] = canonicalTranscript
                ct_field = 'canonicalTranscript_' + jax_id
                canonicalTranscript_graph_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
                editable_statement_writer.writerow([ct_field,canonicalTranscript,edit_date,user_id,canonicalTranscript_graph_id])

                jax_gene_writer.writerow([row[0],canonicalTranscript_graph_id,row[2],row[3],row[4],row[5],EditableSynonymList_graph_id,jax_id])
    return jax_transcript_dict


def migrate_omnigene(load_files_dict,extract_dir, user_id,jax_transcript_dict):
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
    editable_statement_writer = load_files_dict['EditableStatement']['writer']
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
                synonyms = syn_dict[syn_es].split(',')

                EditableSynonymList_graph_id = 'esl_' + now.strftime("%Y%m%d%H%M%S%f")
                for syn in synonyms:
                    synonym_writer.writerow([None,syn,EditableSynonymList_graph_id])
                esyn_field = 'synonyms_' + omnigene_id
                editable_synonym_list_writer.writerow([esyn_field,edit_date,user_id,EditableSynonymList_graph_id])

                canonicalTranscript = ''
                jax_id = row[6]
                if jax_id != None and jax_id in jax_transcript_dict:
                    canonicalTranscript = jax_transcript_dict[jax_id]
                ct_field = 'canonicalTranscript_' + omnigene_id
                canonicalTranscript_graph_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
                editable_statement_writer.writerow([ct_field,canonicalTranscript,edit_date,user_id,canonicalTranscript_graph_id])

                omnigene_writer.writerow([row[0],row[1],row[2],row[3],EditableSynonymList_graph_id,canonicalTranscript_graph_id,row[5],row[6],row[7],omnigene_id])


def get_list_of_files(path: str) -> list:
    json_files = []
    for entry in os.scandir(path):
        if entry.is_file():
            json_files.append(entry.path)
    return json_files


def get_jax_gene_dict(extract_dir)->dict:
    jax_gene_dict = {}
    firstline = True
    with open(extract_dir + 'JaxGene.csv') as csvfile:
        jax_genes = csv.reader(csvfile)
        for row in jax_genes:
            if firstline:
                firstline = False
            else:
                jax_gene_dict[row[4]] = row[6]
    return jax_gene_dict


def get_literature_reference_dict(extract_dir)->dict:
    literature_reference_dict = {}
    firstline = True
    with open(extract_dir + 'LiteratureReference.csv') as csvfile:
        refs = csv.reader(csvfile)
        for row in refs:
            if firstline:
                firstline = False
            else:
                literature_reference_dict[row[0]] = row[10]
    return literature_reference_dict


def get_journal_dict(extract_dir):
    journal_dict = {}
    firstline = True
    with open(extract_dir + 'Journal.csv') as csvfile:
        journals = csv.reader(csvfile)
        for row in journals:
            if firstline:
                firstline = False
            else:
                journal_dict[row[1]] = row[0]

    return journal_dict


def get_author_dict(extract_dir):
    author_dict = {}
    firstline = True
    with open(extract_dir + 'Author.csv') as csvfile:
        authors = csv.reader(csvfile)
        for row in authors:
            if firstline:
                firstline = False
            else:
                author_dict[row[2]] = {'surname':row[0], 'firstInitial':row[1]}

    return author_dict


def preflight_ref(pmid, load_files_dict, data_dict):
    graph_id = None
    ref_by_pmid = data_dict['ref_by_pmid']
    journal_dict = data_dict['journal_dict']
    author_dict = data_dict['author_dict']
    if not pmid in ref_by_pmid:
        reference = get_reference_from_pmid_by_metapub(pmid)
        if reference != None:
            graph_id = 'ref_' + str(pmid)
            journal = reference['journal']
            journal_id = 'journal_' + fix_author_id(journal)
            if not journal_id in journal_dict:
                journal_writer = load_files_dict['Journal']['writer']
                journal_writer.writerow([journal,journal_id])
                journal_dict[journal_id] = journal
            literatureReference_writer = load_files_dict['LiteratureReference']['writer']
            short_ref = ref_name_from_authors_pmid_and_year(reference['authors'], reference['pmid'], reference['year'])
            literatureReference_writer.writerow([reference['pmid'],reference['doi'],reference['title'],journal_id,reference['volume'],reference['first_page'],
                             reference['last_page'],reference['year'],short_ref,reference['abstract'],graph_id])
            author_writer = load_files_dict['Author']['writer']
            literatureReference_author_writer = load_files_dict['LiteratureReference_Author']['writer']
            for author in reference['authors']:
                first, surname = get_authors_names(author)
                author_id = fix_author_id('author_' + surname + '_' + first)
                if not author_id in author_dict:
                    author_writer.writerow([surname, first, author_id])
                    author_dict[author_id] = {'surname':surname, 'firstInitial':first}
                literatureReference_author_writer.writerow([None,author_id,graph_id])
            ref_by_pmid[pmid] = graph_id
    else:
        graph_id = ref_by_pmid[pmid]
    return graph_id


def write_editable_statement(field, editable_statement_writer, loader_id, statement):
    now = datetime.datetime.now()
    es_des_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
    editable_statement_writer.writerow([field, statement, now.strftime("%Y-%m-%d-%H-%M-%S-%f"), loader_id, es_des_id])
    return es_des_id


def read_one_variant_json(path:str)->dict:
    with open(path, 'r') as afile:
        variant_data = json.loads(afile.read())
        ckb_id = str(variant_data['id'])
        full_name = variant_data['fullName']
        variant_type = variant_data['impact']
        protein_effect = variant_data['proteinEffect']
        description = '-'
        if 'geneVariantDescriptions' in variant_data:
            descriptions_array = variant_data['geneVariantDescriptions']
            if len(descriptions_array) > 0:
                description: str = replace_characters(descriptions_array[0]['description'])
                references: list = descriptions_array[0]['references']

        gene_id = variant_data['gene']['id']
        gene = variant_data['gene']['geneSymbol']
        gDot = '-'
        cDot = '-'
        pDot = '-'
        transcript = '-'
        if 'referenceTranscriptCoordinates' in variant_data and variant_data['referenceTranscriptCoordinates'] != None:
            gDot = variant_data['referenceTranscriptCoordinates']['gDna']
            cDot = variant_data['referenceTranscriptCoordinates']['cDna']
            pDot = variant_data['referenceTranscriptCoordinates']['protein']
            transcript = variant_data['referenceTranscriptCoordinates']['transcript']

        variant = {
            'ckb_id':ckb_id,
            'name':full_name,
            'variant_type':variant_type,
            'protein_effect':protein_effect,
            'description':description,
            'references': references,
            'gene_id':str(gene_id),
            'gene':gene,
            'gDot':gDot,
            'cDot':cDot,
            'pDot':pDot,
            'transcript':transcript,
        }
        return variant


def jaxvariant_loader(load_files_dict,data_dict):
    editable_statement_writer = load_files_dict['EditableStatement']['writer']
    editable_protein_effect_writer = load_files_dict['EditableProteinEffect']['writer']
    jaxvariant_writer = load_files_dict['JaxVariant']['writer']
    es_lr_writer = load_files_dict['EditableStatement_LiteratureReference']['writer']

    json_files = get_list_of_files('../data/variants')
    print("num variants=",len(json_files))
    loader_id = data_dict['loader_id']
    jax_gene_dict = data_dict['jax_gene_dict']
    jax_variant_dict = {}
    data_dict['jax_variant_dict'] = jax_variant_dict
    counter = 0
    for json_file in json_files:
        variant: dict = read_one_variant_json(json_file)
        counter += 1
        if (counter % 100 == 0):
            print(counter)
        if variant is not None:
            # print(variant)
            gene_id = variant['gene_id']
            if gene_id in jax_gene_dict:
                graph_id = 'jaxVariant_' + str(variant['ckb_id'])
                jax_variant_dict[str(variant['ckb_id'])] = graph_id
                jaxGene_graph_id = jax_gene_dict[gene_id]

                des_field: str = 'jaxVariantDescription_' + str(variant['ckb_id'])
                es_des_id = write_editable_statement(des_field, editable_statement_writer, loader_id, variant['description'])

                pdot_field: str = 'jaxVariant_pdot_' + str(variant['ckb_id'])
                es_pdot_id = write_editable_statement(pdot_field, editable_statement_writer, loader_id, variant['pDot'])

                cdot_field: str = 'jaxVariant_cdot_' + str(variant['ckb_id'])
                es_cdot_id = write_editable_statement(cdot_field, editable_statement_writer, loader_id, variant['cDot'])

                gdot_field: str = 'jaxVariant_gdot_' + str(variant['ckb_id'])
                es_gdot_id = write_editable_statement(gdot_field, editable_statement_writer, loader_id, variant['gDot'])

                ct_field: str = 'jaxVariant_canonicalTranscript_' + str(variant['ckb_id'])
                es_ct_id = write_editable_statement(ct_field, editable_statement_writer, loader_id, variant['transcript'])

                pe_field: str = 'jaxVariant_protein_effect_' + str(variant['ckb_id'])
                now = datetime.datetime.now()
                es_pe_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
                editable_protein_effect_writer.writerow([pe_field, variant['protein_effect'], now.strftime("%Y-%m-%d-%H-%M-%S-%f"), loader_id, es_pe_id])

                jaxvariant_writer.writerow([variant['name'],es_des_id,variant['ckb_id'],jaxGene_graph_id,es_pdot_id,es_cdot_id,es_gdot_id,es_ct_id,
                                            variant['variant_type'],es_pe_id,graph_id])
                for ref in variant['references']:
                    pmid = ref['pubMedId']
                    if pmid != None:
                        ref_id = preflight_ref(pmid, load_files_dict, data_dict)
                        if ref_id!=None:
                            es_lr_writer.writerow([None, es_des_id, ref_id])


def clinvarvariant_loader(load_files_dict,data_dict):
    print('clinvarvariant_loader')
    print('data_dict:', data_dict)
    data_dict['clinvarvariant_loader'] = 'clinvarvariant_loader'

def hotspotvariant_loader(load_files_dict,data_dict):
    print('hotspotvariant_loader')
    print('data_dict:', data_dict)
    data_dict['hotspotvariant_loader'] = 'hotspotvariant_loader'

def govariant_loader(load_files_dict,data_dict):
    print('govariant_loader')
    print('data_dict:', data_dict)

def variantsnvindel_loader(load_files_dict,data_dict):
    print('variantsnvindel_loader')
    print('data_dict:', data_dict)

def variantregion_loader(load_files_dict,data_dict):
    print('variantregion_loader')
    print('data_dict:', data_dict)

def variantcnv_loader(load_files_dict,data_dict):
    print('variantcnv_loader')
    print('data_dict:',data_dict)

def variantfusion_loader(load_files_dict,data_dict):
    print('variantfusion_loader')
    print('data_dict:', data_dict)

def genomicvariantmarker_loader(load_files_dict,data_dict):
    print('genomicvariantmarker_loader')
    print('data_dict:', data_dict)

def migrate_and_extend(should_migrate,should_generate_sql):
    files_to_copy = ['Author','EditableStatement','EditableStatement_InternetReference',"EditableStatement_LiteratureReference",
                     'InternetReference','Journal','LiteratureReference','LiteratureReference_Author','MyGeneInfoGene',
                     'UniprotEntry','User']

    # files_to_create = ['EditableInt','EditableInt_LiteratureReference','EditableInt_InternetReference',
    #                    'EditableBoolean', 'EditableBoolean_LiteratureReference', 'EditableBoolean_InternetReference',
    #                    'EditableCopyChange', 'EditableCopyChange_LiteratureReference', 'EditableCopyChange_InternetReference',
    #                    'EditableProteinEffect', 'EditableProteinEffect_LiteratureReference', 'EditableProteinEffect_InternetReference',
    #                    'EditableOmniGeneReference', 'EditableOmniGeneReference_LiteratureReference', 'EditableOmniGeneReference_InternetReference',
    #                    'JaxVariant','ClinVarVariant','OncoTreeOccurrence','HotSpotVariant','GOVariant',
    #                    'VariantSNVIndel','VariantRegion','VariantCNV','VariantFusion','GenomicVariantMarker']


    load_dir = get_load_dir()
    extract_dir = get_load_dir() + 'extracted/'
    descriptions_csv_path = '../config/table_descriptions_03_01.csv'
    db_dict = get_schema(descriptions_csv_path)['OmniSeqKnowledgebase2']

    if should_migrate:
        copy_files(files_to_copy, extract_dir, load_dir)

    for file in files_to_copy:
        db_dict.pop(file)

    if not should_migrate:
        db_dict.pop('JaxGene')
        db_dict.pop('OmniGene')
        db_dict.pop('EditableSynonymList')


    load_files_dict = create_load_files_dict(db_dict,load_dir)
    open_file_for_append('EditableStatement', load_dir, load_files_dict)
    open_file_for_append('EditableStatement_LiteratureReference', load_dir, load_files_dict)
    open_file_for_append('EditableStatement_InternetReference', load_dir, load_files_dict)
    open_file_for_append('LiteratureReference', load_dir, load_files_dict)
    open_file_for_append('Author', load_dir, load_files_dict)
    open_file_for_append('LiteratureReference_Author', load_dir, load_files_dict)
    open_file_for_append('Journal', load_dir, load_files_dict)
    open_file_for_append('InternetReference', load_dir, load_files_dict)

    user_id = get_loader_user_id(extract_dir)

    data_dict = {'loader_id':user_id}

    if should_migrate:
        jax_transcript_dict = migrate_jax_gene(load_files_dict, extract_dir, user_id)
        migrate_omnigene(load_files_dict, extract_dir, user_id,jax_transcript_dict)
    else:
        data_dict['jax_gene_dict'] = get_jax_gene_dict(extract_dir)
        data_dict['ref_by_pmid'] = get_literature_reference_dict(extract_dir)
        data_dict['journal_dict'] = get_journal_dict(extract_dir)
        data_dict['author_dict'] = get_author_dict(extract_dir)

    variant_files = ['JaxVariant','ClinVarVariant','HotSpotVariant','GOVariant',
                       'VariantSNVIndel','VariantRegion','VariantCNV','VariantFusion','GenomicVariantMarker']


    for file in variant_files:
        loader_name = file.lower() + '_loader(load_files_dict,data_dict)'
        eval(loader_name )



    close_load_files(load_files_dict)

    if should_generate_sql:
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
    should_migrate = False
    should_generate_sql = False
    migrate_and_extend(should_migrate,should_generate_sql)