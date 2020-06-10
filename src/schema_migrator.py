import csv
import datetime
import os
import shutil
import json
from urllib.parse import urlparse

import mysql.connector
import xml.etree.ElementTree as ET
import xlrd

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

def get_omnigene_dict(extract_dir)->dict:
    omnigene_dict = {}
    firstline = True
    with open(extract_dir + 'OmniGene.csv') as csvfile:
        omni_genes = csv.reader(csvfile)
        for row in omni_genes:
            if firstline:
                firstline = False
            else:
                omnigene_dict[row[0]] = row[9]
    return omnigene_dict


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


def handle_literature_reference_by_pmid(pmid, es_des_id, load_files_dict, data_dict):
    es_lr_writer = load_files_dict['EditableStatement_LiteratureReference']['writer']
    ref_id = preflight_ref(pmid, load_files_dict, data_dict)
    if ref_id != None:
        es_lr_writer.writerow([None, es_des_id, ref_id])


def extract_domain_from_url(url):
    if not url.startswith('http://'):
        url = 'http://' + url
    uri = urlparse(url)
    domain_name = f"{uri.netloc}"
    return domain_name


def handle_internet_reference(web_address, es_id,load_files_dict):
    ir_writer = load_files_dict['InternetReference']['writer']
    es_ir_writer = load_files_dict['EditableStatement_InternetReference']['writer']
    accessedDate = datetime.datetime.now().strftime("%m/%d/%Y")
    shortReference = extract_domain_from_url(web_address) + ' (accessed on:' + accessedDate + ')'
    now = datetime.datetime.now()
    graph_id = 'ref_'  + now.strftime("%Y%m%d%H%M%S%f")
    ir_writer.writerow([accessedDate,web_address,shortReference,graph_id])
    es_ir_writer.writerow([None,es_id,graph_id])


def write_editable_statement(field, editable_statement_writer, loader_id, statement):
    now = datetime.datetime.now()
    es_des_id: str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
    editable_statement_writer.writerow([field, str(statement), now.strftime("%Y-%m-%d-%H-%M-%S-%f"), loader_id, es_des_id])
    return es_des_id


def write_editable_int(field, editable_int_writer, loader_id, the_int):
    now = datetime.datetime.now()
    ei_id: str = 'ei_' + now.strftime("%Y%m%d%H%M%S%f")
    editable_int_writer.writerow([field, int(the_int), now.strftime("%Y-%m-%d-%H-%M-%S-%f"), loader_id, ei_id])
    return ei_id


#  JaxVariant stuff
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
    es_lr_writer = load_files_dict['EditableStatement_LiteratureReference']['writer']
    editable_protein_effect_writer = load_files_dict['EditableProteinEffect']['writer']
    jaxvariant_writer = load_files_dict['JaxVariant']['writer']

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
        # if (counter % 1000 == 0):
        #     break

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
                        handle_literature_reference_by_pmid(pmid, es_des_id, load_files_dict, data_dict)




# ClinVar stuff

def get_cDot(variantName):
    if (len(variantName) > 0 and ':' in variantName):
        vn = variantName.split(':')[1]
        if (' ' in vn):
            variantName = vn.split()[0]
        else:
            variantName = vn
    return variantName


def getSignificanceTuple(sigDict):
    maxVal = 0
    maxName = ''
    explain = ''
    for significance in sigDict:
        if (len(explain) > 0):
            explain += "/"
        explain += significance + "(" + str(sigDict[significance]) + ")"
        if (sigDict[significance] > maxVal):
            maxVal = sigDict[significance]
            maxName = significance
    return maxName, explain

def convert_cdot_to_position(cdot):
    pos = None
    pos_string = ''
    for ch in cdot[2:]:
        if ch.isdigit():
            pos_string = pos_string + ch
        else:
            break
    if len(pos_string) > 0:
        pos = int(pos_string)
    return pos


def getOneVariant(variationArchive):
    clinvar = {
        'variantID': '',
        'gene': '',
        'cDot': '',
        'pDot': '',
        'cdot_pos': 0,
        'pdot_pos': 0,
        'significance': '',
        'signficanceExplanation': '',
    }
    variantName = ''
    sigs = {}
    for simpleAllele in variationArchive.iter('SimpleAllele'):
        if 'VariationID' in simpleAllele.attrib:
            clinvar['variantID'] = simpleAllele.attrib['VariationID']

    for gene in variationArchive.iter('Gene'):
        clinvar['gene'] = gene.attrib['Symbol']
    for name in variationArchive.iter('Name'):
        if (len(variantName) == 0):
            variantName = name.text
    for proteinChange in variationArchive.iter('ProteinChange'):
        if len(clinvar['pDot'])==0:
            clinvar['pDot'] = proteinChange.text
    for clinicalAssertion in variationArchive.iter('ClinicalAssertion'):
        for child in clinicalAssertion:
            if (child.tag == 'Interpretation'):
                for gc in child:
                    if (gc.tag == 'Description'):
                        significance = gc.text.lower()
                        sigs[significance] = sigs.get(significance, 0) + 1
    clinvar['cDot'] = get_cDot(variantName)
    pos = convert_cdot_to_position(clinvar['cDot'])
    if pos is not None:
        clinvar['cdot_pos'] = pos
        clinvar['pdot_pos'] = int((int(pos) + 2) / 3)

    clinvar['significance'], clinvar['signficanceExplanation'] = getSignificanceTuple(sigs)

    return clinvar

def clinvarvariant_loader(load_files_dict,data_dict):
    print('clinvarvariant_loader')
    editable_statement_writer = load_files_dict['EditableStatement']['writer']
    clinvarvariant_writer = load_files_dict['ClinVarVariant']['writer']
    loader_id = data_dict['loader_id']
    counter = 0
    id_dict = {}
    clinvar_variant_dict = {}
    data_dict['clinvar_variant_dict'] = clinvar_variant_dict
    for event, elem in ET.iterparse('../data/ClinVarVariationRelease_2020-05.xml'):
        if event == 'end':
            if elem.tag == 'VariationArchive':
                clinvar = getOneVariant(elem)
                id = str(clinvar['variantID'])
                if id in id_dict:
                    now = datetime.datetime.now()
                    id = id +  now.strftime("%Y%m%d%H%M%S%f")
                id_dict[id] = id
                graphql = 'clinvar_variant_' + id
                if clinvar['pDot'] != '':
                    clinvar_variant_dict[clinvar['pDot']] = graphql

                gene_field: str = 'clinvar_variant_gene_' + id
                es_gene_id = write_editable_statement(gene_field, editable_statement_writer, loader_id, clinvar['gene'])

                pdot_field: str = 'clinvar_variant_pdot_' + id
                es_pdot_id = write_editable_statement(pdot_field, editable_statement_writer, loader_id, clinvar['pDot'])

                cdot_field: str = 'clinvar_variant_cdot_' + id
                es_cdot_id = write_editable_statement(cdot_field, editable_statement_writer, loader_id, clinvar['cDot'])

                gdot_field: str = 'clinvar_variant_significance_' + id
                es_significance_id = write_editable_statement(gdot_field, editable_statement_writer, loader_id, clinvar['significance'])

                ct_field: str = 'clinvar_variant_signficanceExplanation_' + id
                es_signficanceExplanation_id = write_editable_statement(ct_field, editable_statement_writer, loader_id, clinvar['signficanceExplanation'])

                if not 'HAPLOTYPE' in clinvar['cDot'] and len(clinvar['cDot'])<100:
                    clinvarvariant_writer.writerow([clinvar['variantID'],es_gene_id,es_pdot_id,es_cdot_id,es_significance_id,es_signficanceExplanation_id,graphql])
                counter += 1
                if (counter % 1000 == 0):
                    print(counter)
                # if (counter % 10000 == 0):
                #     break
                elem.clear()  # discard the element

# Hot spot stuff
def handle_occurrences(occurrences, onco_dict):
    occurrences_list = []
    occurrence_array = occurrences.split('|')
    for item in occurrence_array:
        vals = item.split(':')
        disease = vals[0]
        total = int(vals[1])
        has_variant = int(vals[2])
        ratio = 100.0 * (has_variant / total)
        onco_tree_occurrence = {'percentOccurrence': '%.1f%%' % ratio,
                                'perThousand': round(10.0 * ratio),
                                'occurences': has_variant,
                                'totalSamples': total}
        if disease in onco_dict:
            onco_tree_occurrence['disease'] = onco_dict[disease]['name']
            onco_tree_occurrence['oncoTreeCode'] = onco_dict[disease]['code']
        else:
            onco_tree_occurrence['disease'] = disease
            onco_tree_occurrence['oncoTreeCode'] = disease
        occurrences_list.append(onco_tree_occurrence)
    return occurrences_list

def read_snv_hotspot(onco_dict):
    hot_spots = []
    with open('../data/hotspots_v2/SNV-hotspots-Table 1.tsv') as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            # print(row)
            gene = row['Hugo_Symbol']
            referenceAminoAcid = row['Reference_Amino_Acid'].split(':')[0]
            variantAminoAcid = row['Variant_Amino_Acid'].split(':')[0]
            begin = row['Amino_Acid_Position']
            if not begin.isnumeric():
                position = 0
            else:
                position = int(begin)
            if referenceAminoAcid == 'splice':
                name = gene + ' ' + begin
            else:
                name = gene +  ' ' + referenceAminoAcid + begin + variantAminoAcid

            occurrences = handle_occurrences(row['Detailed_Cancer_Types'],onco_dict)
            hot_spot = {'name':name, 'gene':gene,
                        'referenceAminoAcid':referenceAminoAcid,
                        'variantAminoAcid':variantAminoAcid,
                        'begin':begin,
                        'end':begin,
                        'position':position,
                        'occurrences':occurrences}
            hot_spots.append(hot_spot)
    return hot_spots

def read_indel_hotspot(onco_dict):
    hot_spots = []
    with open('../data/hotspots_v2/INDEL-hotspots-Table 1.tsv') as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            # print(row)
            gene = row['Hugo_Symbol']

            variant = row['Variant_Amino_Acid'].split(':')[0]
            name = gene + ' ' + variant
            referenceAminoAcid = variant[:1]
            variantAminoAcid = ''
            if '_' in variant:
                variantAminoAcid = variant.split('_')[1][:1]
            if '-' in row['Amino_Acid_Position']:
                pos = row['Amino_Acid_Position'].split('-')
                begin = pos[0]
                end = pos[1]
            else:
                begin = row['Amino_Acid_Position']
                end = begin
            position = int(begin)
            if position < 0:
                position = 0

            occurrences = handle_occurrences(row['Detailed_Cancer_Types'],onco_dict)
            hot_spot = {'name':name, 'gene':gene,
                        'referenceAminoAcid':referenceAminoAcid,
                        'variantAminoAcid':variantAminoAcid,
                        'begin':begin,
                        'end':end,
                        'position':position,
                        'occurrences':occurrences}
            hot_spots.append(hot_spot)
    return hot_spots

def get_oncotree_dict():
    onco_dict = {}
    unknown = {'code':'unk', 'name':'Unknown'}
    onco_dict['unk'] = unknown
    f = open('../data/onctoree.json', "r")
    data = json.loads(f.read())
    f.close()
    for item in data:
        onco_dict[item['code'].lower()] = item
        for h in item['history']:
            onco_dict[h.lower()] = item
        for p in item['precursors']:
            onco_dict[p.lower()] = item
    return onco_dict


def hotspotvariant_loader(load_files_dict,data_dict):
    print('hotspotvariant_loader')
    editable_statement_writer = load_files_dict['EditableStatement']['writer']
    editable_int_writer = load_files_dict['EditableInt']['writer']
    hotspotvariant_writer = load_files_dict['HotSpotVariant']['writer']
    oncotreeoccurrence_writer = load_files_dict['OncoTreeOccurrence']['writer']
    loader_id = data_dict['loader_id']

    hotspot_variant_dict = {}
    data_dict['hotspot_variant_dict'] = hotspot_variant_dict


    onco_dict = get_oncotree_dict()
    hot_spots = read_snv_hotspot(onco_dict)
    hot_spots.extend(read_indel_hotspot(onco_dict))
    for hot_spot in hot_spots:
        id = hot_spot['name'].replace(' ', '_').replace('-', '_').replace('*', '')
        graph_id = 'hot_spot_' + id
        hotspot_variant_dict[hot_spot['name']] = graph_id
        begin_field: str = 'hot_spot_begin_' + id
        es_begin_id = write_editable_statement(begin_field, editable_statement_writer, loader_id, hot_spot['begin'])

        end_field: str = 'hot_spot_end_' + id
        es_end_id = write_editable_statement(end_field, editable_statement_writer, loader_id, hot_spot['end'])

        position_field: str = 'hotspotvariant_position_' + id
        ei_position_id: str = write_editable_int(position_field, editable_int_writer, loader_id, int(hot_spot['position']))

        # hot_spot['name'],hot_spot['gene'],hot_spot['referenceAminoAcid'],hot_spot['variantAminoAcid'],
        #                                                    hot_spot['begin'],hot_spot['end'],str(hot_spot['position']),graph_id

        hotspotvariant_writer.writerow([hot_spot['name'],hot_spot['gene'],hot_spot['referenceAminoAcid'],hot_spot['variantAminoAcid'],
                                        es_begin_id,es_end_id,ei_position_id,graph_id])
        # hot_spot_occurrence['disease'],hot_spot_occurrence['oncoTreeCode'],hot_spot_occurrence['percentOccurrence'],
        #                                                    str(hot_spot_occurrence['perThousand']),str(hot_spot_occurrence['occurences']),str(hot_spot_occurrence['totalSamples']),hot_spot_id
        for occurrence in hot_spot['occurrences']:
            occurrence_graph_id = 'oncotree_occurrence_' + graph_id + '_' + occurrence['oncoTreeCode']
            disease_graph_id = 'oncotree_disease_' + occurrence['oncoTreeCode']
            percentOccurrence_field: str = 'percentOccurrence_' + disease_graph_id
            es_percentOccurrence_id = write_editable_statement(percentOccurrence_field, editable_statement_writer, loader_id, occurrence['percentOccurrence'])
            occurrences_field: str = 'occurrences_' + disease_graph_id
            ei_occurrences_id: str = write_editable_int(occurrences_field, editable_int_writer, loader_id, int(occurrence['occurences']))
            totalSamples_field: str = 'totalSamples_' + disease_graph_id
            ei_totalSamples_id: str = write_editable_int(totalSamples_field, editable_int_writer, loader_id, int(occurrence['totalSamples']))

            oncotreeoccurrence_writer.writerow([disease_graph_id,occurrence['oncoTreeCode'],es_percentOccurrence_id,ei_occurrences_id,ei_totalSamples_id,
                                                str(occurrence['perThousand']),graph_id,occurrence_graph_id])

#  go variant stuff
def read_one_go_json(path:str)->dict:
    with open(path, 'r') as afile:
        go_data = json.loads(afile.read())
        results = go_data['results']
        return results

def govariant_loader(load_files_dict,data_dict):
    print('govariant_loader')
    json_files = get_list_of_files('../data/GO_alterations')
    print("num files=", len(json_files))
    loader_id = data_dict['loader_id']
    editable_statement_writer = load_files_dict['EditableStatement']['writer']
    govariant_writer = load_files_dict['GOVariant']['writer']
    jax_dict = data_dict['jax_variant_dict']
    go_variant_dict = {}
    data_dict['go_variant_dict'] = go_variant_dict

    counter = 0
    go_dict = {}
    for json_file in json_files:
        print(json_file)
        alteration_array = read_one_go_json(json_file)
        for alteration in alteration_array:
            if 'gene' in alteration:
                go_id = alteration['id']
                if not go_id in go_dict:
                    if 'mutation_type' in alteration:
                        mutation_type = alteration['mutation_type']
                    else:
                        mutation_type = ''
                    graph_id = 'go_variant_' + go_id
                    variant_name = alteration['name']
                    go_variant_dict[variant_name] = graph_id
                    gene_name = alteration['gene']
                    jax_variant_id = None
                    if 'codes' in alteration:
                        for code in alteration['codes']:
                            if code.startswith('JAX'):
                                jax_id = code[12:]
                                if jax_id in jax_dict:
                                    jax_variant_id = jax_dict[jax_id]
                                # variant_name,gene_name,go_id,mutationType,jax_variant_id,graph_id
                                break
                    name_field: str = 'go_variant_name_' + go_id
                    es_name_id = write_editable_statement(name_field, editable_statement_writer, loader_id, variant_name)
                    gene_field: str = 'go_variant_gene_' + go_id
                    es_gene_id = write_editable_statement(gene_field, editable_statement_writer, loader_id, gene_name)
                    mutationType_field: str = 'go_variant_mutationType_' + go_id
                    es_mutationType_id = write_editable_statement(mutationType_field, editable_statement_writer, loader_id, mutation_type)

                    govariant_writer.writerow([es_name_id,es_gene_id,go_id,es_mutationType_id,jax_variant_id,graph_id])
                    go_dict[go_id] = graph_id
                    counter += 1
                    if (counter % 100 == 0):
                        print(counter)


def create_variants_data_dicts(data_dict):
    protein_effect_dict = {
        'gain of function': 'GainOfFunction', 'gain of function - predicted': 'GainOfFunctionPredicted',
        'loss of function': 'LossOfFunction', 'loss of function - predicted': 'LossOfFunctionPredicted',
        'no effect': 'NoEffect', 'unknown': 'Unknown'}
    data_dict['protein_effect_dict'] = protein_effect_dict
    variants_data_dict = {}
    data_dict['variants_data_dict'] = variants_data_dict
    wb = xlrd.open_workbook(filename='data/OCP_Lev1_2_Variants_w_exon.xlsx')
    sheet = wb.sheet_by_index(0)
    for row_idx in range(1, sheet.nrows):
        geneName = sheet.cell_value(row_idx, 0)
        variantName = sheet.cell_value(row_idx, 1)
        copyChange = sheet.cell_value(row_idx, 3)
        fusionDesc = sheet.cell_value(row_idx, 4)
        regionType = sheet.cell_value(row_idx, 9)
        regionNumber = sheet.cell_value(row_idx, 12)
        if isinstance(regionNumber,float):
            regionNumber = str(int(regionNumber))
        regionAAChange = sheet.cell_value(row_idx, 13)
        regionIndelType= sheet.cell_value(row_idx, 14)
        exon = sheet.cell_value(row_idx, 18)
        exon_number = 0;
        if isinstance(exon,float):
            exon_number = int(exon)
            exon = str(int(exon))
        pdot = ''
        exon3 = ''
        exon5 = ''
        if fusionDesc != '':
            variantType = 'Fusion'
            if exon != '':
                pdot = 'Exon ' + exon + ' Deletion'
                exon5 = str(exon_number - 1)
                exon3 = str(exon_number + 1)
            elif 'rearrange' in variantName:
                pdot = 'rearrangement'
            else:
                pdot = 'fusion'
        elif copyChange != '':
            variantType = 'CNV'
            pdot = copyChange
        elif regionIndelType != '':
            if regionAAChange=='':
                variantType = 'Region'
                pdot = regionType + ' ' + regionNumber + ' ' + regionIndelType
            else:
                variantType = 'Indel'
                pdot = regionAAChange
        else:
            if regionAAChange=='':
                variantType = 'Region'
                if regionType == 'Gene':
                    pdot = ' mutation'
                else:
                    pdot = regionType + ' ' + regionNumber + ' mutation'
            else:
                variantType = 'SNV'
                pdot = regionAAChange
        gene_plus_pdot = geneName + ' ' + pdot
        variant = {'geneName':geneName,'variantName':variantName,'variantType':variantType, 'copyChange':copyChange, 'fusionDesc':fusionDesc,
                   'regionType':regionType, 'regionNumber':regionNumber,'pdot':pdot,'exon5':exon5,'exon3':exon3, 'gene_plus_pdot':gene_plus_pdot }
        variants_data_dict[gene_plus_pdot] = variant



def variantsnvindel_loader(load_files_dict,data_dict):
    print('variantsnvindel_loader')
    if not 'variants_data_dict' in data_dict:
        create_variants_data_dicts(data_dict)
    variants_data_dict = data_dict['variants_data_dict']

def variantregion_loader(load_files_dict,data_dict):
    print('variantregion_loader')
    if not 'variants_data_dict' in data_dict:
        create_variants_data_dicts(data_dict)
    variants_data_dict = data_dict['variants_data_dict']

def variantcnv_loader(load_files_dict,data_dict):
    print('variantcnv_loader')
    if not 'variants_data_dict' in data_dict:
        create_variants_data_dicts(data_dict)
    variants_data_dict = data_dict['variants_data_dict']

def variantfusion_loader(load_files_dict,data_dict):
    print('variantfusion_loader')
    if not 'variants_data_dict' in data_dict:
        create_variants_data_dicts(data_dict)
    variants_data_dict = data_dict['variants_data_dict']

def genomicvariantmarker_loader(load_files_dict,data_dict):
    print('genomicvariantmarker_loader')
    if not 'variants_data_dict' in data_dict:
        create_variants_data_dicts(data_dict)
    variants_data_dict = data_dict['variants_data_dict']


def migrate_and_extend(should_migrate,should_generate_sql):
    print(datetime.datetime.now().strftime("%H:%M:%S"))
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

    data_dict['jax_gene_dict'] = get_jax_gene_dict(extract_dir)
    data_dict['ref_by_pmid'] = get_literature_reference_dict(extract_dir)
    data_dict['journal_dict'] = get_journal_dict(extract_dir)
    data_dict['author_dict'] = get_author_dict(extract_dir)
    data_dict['omnigene_dict'] = get_omnigene_dict(extract_dir)

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
    print(datetime.datetime.now().strftime("%H:%M:%S"))


if __name__ == "__main__":
    should_migrate = False
    should_generate_sql = False
    migrate_and_extend(should_migrate,should_generate_sql)