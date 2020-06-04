import sys
import os
import time
import datetime
import requests
os.environ['NCBI_API_KEY'] = 'cde5c1a63fa16711994bfe74b858747cbb08'
from metapub import PubMedFetcher
import re
from neo4j import GraphDatabase
import unicodedata as ud

cache = '/Users/mglynias/Documents/GitHub/OmniSeqKnowledgebase_populate/cache'


def rmdiacritics(char):
    '''
    Return the base character of char, by "removing" any
    diacritics like accents or curls and strokes and the like.  ð
    '''
    if char == "æ":
        return 'ae'
    elif char == 'ß':
        return 's'
    elif char == 'ð':
        return 'd'
    else:
        desc = ud.name(char)
        cutoff = desc.find(' WITH ')
        if cutoff != -1:
            desc = desc[:cutoff]
            try:
                char = ud.lookup(desc)
            except KeyError:
                pass  # removing "WITH ..." produced an invalid name
        return char

def remove_accents(input_str):
    nfkd_form = ud.normalize('NFKD', input_str)
    return u"".join([rmdiacritics(c) for c in input_str])


def send_query(query:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    responseBody = ''
    try:
        response = requests.request("POST", url, headers=headers, json={'query': query})
        if not response.ok:
            response.raise_for_status()
            os.system("say another error")

            sys.exit()

        responseBody: str = response.json()
        # print(responseBody)
        if 'errors' in responseBody:
            print(responseBody)
            os.system("say another error")
            sys.exit()
    except requests.exceptions.RequestException as err:
        print('error in request')
        print(err)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    return responseBody


def send_schema_request(query:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    responseBody = ''
    try:
        response = requests.request("POST", url, headers=headers, json={'query': query})
        if not response.ok:
            response.raise_for_status()
            os.system("say another error")

            sys.exit()

        responseBody: str = response.json()
        # print(responseBody)
        if 'errors' in responseBody:
            print(responseBody)
            os.system("say another error")
            sys.exit()
    except requests.exceptions.RequestException as err:
        print('error in request')
        print(err)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    return responseBody


def send_mutation(mutation_payload:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    mutation_payload = '{"query":"mutation {' + mutation_payload + '}"}'
    # print(mutation_payload)
    responseBody = ''
    try:
        response = requests.request("POST", url, headers=headers, data = mutation_payload.encode('utf-8'))
        if not response.ok:
            response.raise_for_status()
            print(mutation_payload)
            os.system("say another error")

            sys.exit()

        responseBody: str = response.json()
        # print(responseBody)
        if 'errors' in responseBody:
            print(mutation_payload)
            print(responseBody)
            os.system("say another error")
            sys.exit()
    except requests.exceptions.RequestException as err:
        print('error in request')
        print(err)
        print(mutation_payload)
        print(response.text)
        os.system("say another error")
        sys.exit()
    except UnicodeEncodeError as err:
        print('UnicodeEncodeError')
        print(err)
        print(mutation_payload)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    # print(responseBody)
    return responseBody


def get_editor_id(name:str,server:str)->str:
    id = None
    query = f'{{ User(name:"{name}"){{id}} }}'

    response = send_query(query, server)
    if len(response['data']['User'])>0:
        result = response['data']['User'][0]
        id = result['id']
    return id


def get_editor_ids(server:str)->dict:
    id = None
    query = f'{{ User{{id,name}} }}'
    user_dict = {}
    response = send_query(query, server)
    if len(response['data']['User'])>0:

        for u in response['data']['User']:
            name = u['name']
            id = u['id']
            user_dict[name] = id
    return user_dict

def get_jax_descriptions(server)->dict:
    jax_dict: dict = {}
    query = f'{{ JaxGene  {{ id, name, description {{ statement,field }} }} }}'
    response = send_query(query, server)
    if len(response['data']['JaxGene'])>0:
        for item in response['data']['JaxGene']:
            id: str = item['id']
            name: str = item['name']
            description: str = item['description']
            statement: str = description['statement']
            field: str = description['field']
            jax_dict[name] = {'id':id, 'statement':statement, 'field':field}
    return jax_dict


def get_jax_gene_ids(server)->dict:
    jax_dict: dict = {}
    query = f'{{ JaxGene  {{ id, name }} }}'
    response = send_query(query, server)
    if len(response['data']['JaxGene'])>0:
        for item in response['data']['JaxGene']:
            id: str = item['id']
            name: str = item['name']
            jax_dict[name] = id
    return jax_dict


def get_dict_from_omnigene_es_fragment(es):
    id: str =es['id']
    statement: str = es['statement']
    field: str = es['field']
    es_dict = {'id': id, 'statement': statement, 'field': field}
    return es_dict


def get_omnigene_descriptions(server)->dict:
    omnigene_dict: dict = {}
    query = f'{{ OmniGene  {{ id, name, geneDescription {{ id, statement,field }}, oncogenicCategory {{ id, statement,field }}, synonymsString {{ id, statement,field }} }} }}'
    response = send_query(query, server)
    if len(response['data']['OmniGene'])>0:
        for item in response['data']['OmniGene']:
            id: str = item['id']
            name: str = item['name']
            description = get_dict_from_omnigene_es_fragment(item['geneDescription'])
            oncogenic_category = get_dict_from_omnigene_es_fragment(item['oncogenicCategory'])
            synonyms = get_dict_from_omnigene_es_fragment(item['synonymsString'])
            omnigene_dict[name] = {'id':id, 'description':description, 'oncogenic_category':oncogenic_category, 'synonyms':synonyms}
    return omnigene_dict


def get_authors(server:str)->dict:
    author_dict: dict = {}
    query = f'{{ Author  {{ id,surname,first_initial }} }}'
    response = send_query(query, server)
    if len(response['data']['Author'])>0:
        for item in response['data']['Author']:
            id = item['id']
            surname = item['surname']
            first_initial = item['first_initial']
            key = fix_author_id(surname + '_' + first_initial)
            author_dict[key] = id
    return author_dict


def get_literature_references(server:str)->dict:
    reference_dict: dict = {}
    query = f'{{ LiteratureReference  {{ id,PMID }} }}'
    response = send_query(query, server)
    if len(response['data']['LiteratureReference'])>0:
        for item in response['data']['LiteratureReference']:
            id = item['id']
            pmid = item['PMID']
            reference_dict[pmid] = id
    return reference_dict


def get_journals(server:str)->dict:
    journal_dict: dict = {}
    query = f'{{ Journal  {{ id,name }} }}'
    response = send_query(query, server)
    if len(response['data']['Journal'])>0:
        for item in response['data']['Journal']:
            id = item['id']
            name = item['name']
            journal_dict[name] = id
    return journal_dict


def get_omnigenes(server:str)->dict:
    omnigene_dict: dict = {}
    query = f'{{ OmniGene  {{ id,name }} }}'
    response = send_query(query, server)
    if len(response['data']['OmniGene'])>0:
        for item in response['data']['OmniGene']:
            id = item['id']
            name = item['name']
            omnigene_dict[name] = id
    return omnigene_dict


def replace_characters(a_string: str):
    if a_string is not None:
        a_string = a_string.replace('α', 'alpha')
        a_string = a_string.replace(u"\u0391", 'A')
        a_string = a_string.replace('β', 'beta')
        a_string = a_string.replace('γ', 'gamma')
        a_string = a_string.replace('δ', 'delta')
        a_string = a_string.replace('ε', 'epsilon')
        a_string = a_string.replace('ζ', 'zeta')
        a_string = a_string.replace('η', 'eta')
        a_string = a_string.replace('θ', 'theta')
        a_string = a_string.replace('ι', 'iota')
        a_string = a_string.replace('ɩ', 'iota')
        a_string = a_string.replace('κ', 'kappa')
        a_string = a_string.replace('λ', 'lamda')
        a_string = a_string.replace('μ', 'mu')
        a_string = a_string.replace('ν', 'nu')
        a_string = a_string.replace('π', 'pi')
        a_string = a_string.replace('ρ', 'rho')
        a_string = a_string.replace('σ', 'sigma')
        a_string = a_string.replace('χ', 'chi')

        a_string = a_string.replace('ω', 'omega')
        a_string = a_string.replace(u"\u0394", 'delta')

        a_string = a_string.replace(u"\u03c5", 'upsilon')
        a_string = a_string.replace(u"\u03a5", 'Upsilon')
        a_string = a_string.replace('Ψ', 'Psi')
        a_string = a_string.replace('Ω', 'Omega')

        a_string = a_string.replace(u"\u025b", 'e')
        a_string = a_string.replace(u"\u0190", 'e')
        a_string = a_string.replace(u"\u223c", '~')
        a_string = a_string.replace(u"\u301c", '~')


        a_string = a_string.replace("á", "a")
        a_string = a_string.replace("à", "a")
        a_string = a_string.replace("ä", "a")
        a_string = a_string.replace("å", "a")
        a_string = a_string.replace("ã", "a")
        a_string = a_string.replace("â", "a")
        a_string = a_string.replace("ą", "a")
        a_string = a_string.replace("æ", "ae")

        a_string = a_string.replace("ç", "c")
        a_string = a_string.replace("č", "c")
        a_string = a_string.replace("ć", 'c')
        #
        a_string = a_string.replace("ě", "e")
        a_string = a_string.replace("ė", "e")
        a_string = a_string.replace("ę", "e")
        a_string = a_string.replace("é", "e")
        a_string = a_string.replace("è", "e")
        a_string = a_string.replace("ë", "e")
        a_string = a_string.replace("ê", "e")
        #
        a_string = a_string.replace("ﬁ", "fi")
        a_string = a_string.replace("ğ", "g")

        a_string = a_string.replace("í", "i")
        a_string = a_string.replace("ì", "i")
        a_string = a_string.replace("î", "i")
        a_string = a_string.replace("ï", "i")

        a_string = a_string.replace("ń", "n")
        a_string = a_string.replace("ň", "n")
        a_string = a_string.replace("ñ", "n")

        a_string = a_string.replace("ő", "o")
        a_string = a_string.replace("õ", "o")
        a_string = a_string.replace("ö", "o")
        a_string = a_string.replace("ó", "o")
        a_string = a_string.replace("ò", "o")
        a_string = a_string.replace("ô", "o")
        a_string = a_string.replace("ø", "o")

        a_string = a_string.replace("ř", "r")

        a_string = a_string.replace("ş", "s")
        a_string = a_string.replace("ś", "s")
        a_string = a_string.replace("š", "s")
        a_string = a_string.replace("Š", "S")
        a_string = a_string.replace("Ş", "S")
        a_string = a_string.replace("ß", "s")

        a_string = a_string.replace("ť", "t")
        a_string = a_string.replace("ů", "u")
        a_string = a_string.replace("ü", "u")
        a_string = a_string.replace("ū", "u")
        a_string = a_string.replace("ù", "u")
        a_string = a_string.replace("ú", "u")

        a_string = a_string.replace("ÿ", "y")
        a_string = a_string.replace("ý", "y")
        a_string = a_string.replace("ż", "z")
        a_string = a_string.replace("ź", "z")
        a_string = a_string.replace("ž", "z")

        a_string = a_string.replace("’", "")
        a_string = a_string.replace('"', '')
        a_string = a_string.replace('\\', ' ')
        a_string = a_string.replace(u"\u2216", ' ')

        a_string = a_string.replace(u"\u201c", '')
        a_string = a_string.replace(u"\u201d", '')
        a_string = a_string.replace(u"\u2018", '')
        a_string = a_string.replace(u"\u2019", '')
        a_string = a_string.replace(u"\u05f3", '')
        a_string = a_string.replace(u"\u2032", '_')
        a_string = a_string.replace(u"\u2020", '*')
        a_string = a_string.replace(u"\u0142", '')
        a_string = a_string.replace(u"\u202f", ' ')
        a_string = a_string.replace(u"\u200a", ' ')
        a_string = a_string.replace(u"\u2002", ' ')
        a_string = a_string.replace('→', '->')
        a_string = a_string.replace(u"\u2012", '-')
        a_string = a_string.replace(u"\u207b", '-')
        a_string = a_string.replace(u"\uff0c", ', ')

        a_string = a_string.replace(u"\u207a", '+')
        a_string = a_string.replace(u"\u2011", '-')
        a_string = a_string.replace(u"\u2013", '-')
        a_string = a_string.replace(u"\u2014", '-')
        a_string = a_string.replace(u"\u2044", '/')
        a_string = a_string.replace(u"\u2122", 'TM')
        a_string = a_string.replace(u"\u2005", ' ')
        a_string = a_string.replace(u"\u2009", ' ')
        a_string = a_string.replace(u"\u0131", 'i')
        a_string = a_string.replace(u"\u2081", '1')
        a_string = a_string.replace(u"\u2082", '2')
        a_string = a_string.replace(u"\u2265", '>=')
        a_string = a_string.replace(u"\u2264", '<=')
        a_string = a_string.replace(u"\u2264", '<=')
        a_string = a_string.replace(u"\u226b", ' >>')
        a_string = a_string.replace(u"\u2248", ' =')
        a_string = a_string.replace('\t',' ')
        a_string = a_string.replace('\r','')
        a_string = a_string.replace('\n','')
        a_string = a_string.replace('⁸⁸','')
        a_string = a_string.replace('⁹⁰','')
        a_string = a_string.replace('Ⅱ','II')
        a_string = a_string.replace('Ⅰ','I')
        a_string = a_string.replace('&', '')
    return a_string

def fix_author_id(id:str)->str:
    id = id.lower()
    id = remove_accents(id)
    id = id.replace(" ", "")
    id = id.replace(":", "")
    id = id.replace(",", "")
    id = id.replace("(", "")
    id = id.replace(")", "")
    id = id.replace("<sup>®<_sup>","")
    id = id.replace("<", "")
    id = id.replace(">", "")
    id = id.replace("®", "_")
    id = id.replace("-", "_")
    id = id.replace("'", "_")
    id = id.replace("ʼ", "_")
    id = id.replace("ʼ", "_")
    id = id.replace(".", "_")
    id = id.replace("/", "_")

    return id

def PMID_extractor(text:str)->list:
    pattern = r'PMID:\s+\d{8}'
    matches = re.findall(pattern,text)
    pmids = []
    for match in matches:
        if match not in pmids:
            pmids.append(match)
    return pmids

def PubMed_extractor(text:str)->list:
    pattern = r'PubMed:\d{8}'
    matches = re.findall(pattern,text)
    pmids = []
    for match in matches:
        match = match[7:]
        if match not in pmids:
            pmids.append(match)
    return pmids

def get_reference_from_pmid_by_metapub(pmid:str)->dict:
    fetch = PubMedFetcher(cachedir=cache)
    reference = None
    try:
        time.sleep(0.34)
        article = fetch.article_by_pmid(pmid)
        reference = {'journal':article.journal,
                     'authors': article.authors,
                     'issue':article.issue,
                     'first_page':article.first_page,
                     'last_page': article.last_page,
                     'volume':article.volume,
                     'year': str(article.year),
                     'abstract': replace_characters(article.abstract),
                     'title': replace_characters(article.title),
                     'doi': article.doi,
                     'pmid': article.pmid
                     }
    except:
        print('*** Bad PMID:',pmid)

    return reference

def get_authors_names(author):
    l = author.split()
    surname = replace_characters(l[0])
    first = '-'
    if len(l)>1:
        first = replace_characters(l[1])
    return first, surname

def ref_name_from_authors_pmid_and_year(authors, pmid, year):
    s = ''
    if len(authors)>0:
        first, surname = get_authors_names(authors[0])
        if len(authors) == 1:
            s += surname + ' ' + year
        elif len(authors) == 2:
            first2, surname2 = get_authors_names(authors[1])
            s += surname + ' & '+ surname2 + ' ' + year
        else:
            s += surname + ' et al. ' + year
    else:
        s += 'no_authors ' + year
    s += ' (PMID:' + pmid + ')'
    return s


def create_reference_mutation(ref_id, ref):
    ref_name = ref_name_from_authors_pmid_and_year(ref['authors'], ref['pmid'], ref['year'])
    s = f'''{ref_id}: createLiteratureReference(id: \\"{ref_id}\\", abstract: \\"{ref['abstract']}\\", shortReference: \\"{ref_name}\\", title: \\"{ref['title']}\\", volume: \\"{ref['volume']}\\", firstPage: \\"{ref['first_page']}\\", lastPage: \\"{ref['last_page']}\\", publicationYear: \\"{ref['year']}\\", DOI: \\"{ref['doi']}\\", PMID: \\"{ref['pmid']}\\"),'''
    return s


def create_author_mutation(id,surname,first):
    s = f'''{id}: createAuthor(firstInitial: \\"{first}\\" , id: \\"{id}\\",surname: \\"{surname}\\"),'''
    return s


def create_journal_mutation(journal, journal_id):
    s = f'''{journal_id}: createJournal(id: \\"{journal_id}\\",name: \\"{journal}\\"),'''
    return s


def create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id):
    id = ref_id + '_' + journal_id
    s = f'{id}: addLiteratureReferenceJournal(id:\\"{ref_id}\\", journal:\\"{journal_id}\\"),'
    return s

def create_AddLiteratureReferenceAuthors_mutation(ref_id, authors):
    id = 'author_' +ref_id
    author_string = '['
    for a in authors:
        if len(author_string)>1:
            author_string += ","
        author_string += '\\"' + a + '\\"'
    author_string += ']'
    s = f'{id}: addLiteratureReferenceAuthors(id:\\"{ref_id}\\", authors:{author_string}),'

    return s

def write_references(es_id:str, description:str,pmid_extractor:callable,reference_dict:dict,journal_dict:dict,author_dict:dict)->str:
    s: str = ''
    references = []
    pmids = pmid_extractor(description)
    for pmid in pmids:
        reference = get_reference_from_pmid_by_metapub(pmid)
        if reference != None:
            references.append(reference)
    if len(references)>0:
        reference_string = '['
        for r in references:
            pubmed = r['pmid']

            if pubmed not in reference_dict:
                ref_id = 'ref_' + pubmed
                s += create_reference_mutation(ref_id, r)
                reference_dict[pubmed] = ref_id
                journal = r['journal']
                if journal not in journal_dict:
                    journal_id = 'journal_' + fix_author_id(journal)
                    s += create_journal_mutation(journal, journal_id)
                    journal_dict[journal] = journal_id
                else:
                    journal_id = journal_dict[journal]
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in r['authors']:
                    first, surname = get_authors_names(author)
                    key = fix_author_id(surname + '_' + first)
                    if key not in author_dict:
                        author_id = 'author_' + surname + '_' + first
                        author_id = fix_author_id(author_id)
                        s += create_author_mutation(author_id, surname, first)
                        author_dict[key] = author_id
                    else:
                        author_id = author_dict[key]
                    authors.append(author_id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)
            else:
                ref_id = reference_dict[pubmed]
            reference_string += '\\"' + ref_id + '\\",'
        reference_string += ']'
        ref_id2 = 'esref_' + datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        s += f'{ref_id2}: addEditableStatementReferences(id:\\"{es_id}\\", references:{reference_string}),'
    return s

def createEditableStatement(statement:str, field:str, editor_id:str,pmid_extractor:callable,reference_dict:dict,journal_dict:dict,author_dict:dict) -> (str,str):
    now = datetime.datetime.now()
    edit_date:str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
    id:str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
    ede_id:str = 'ese_' + now.strftime("%Y%m%d%H%M%S%f")
    s = f'''{id} : createEditableStatement(deleted: false, edit_date: \\"{edit_date}\\", field: \\"{field}\\", id: \\"{id}\\",statement: \\"{statement}\\"),'''
    s += f'{ede_id}: addEditableStatementEditor(editor:[\\"{editor_id}\\"], id:\\"{id}\\" ),'
    s += write_references(id,statement,pmid_extractor,reference_dict,journal_dict,author_dict)
    return s, id

def createEditableStatement_with_date(statement:str, field:str, editor_id:str,edit_date:str,pmid_extractor:callable,reference_dict:dict,journal_dict:dict,author_dict:dict) -> (str,str):
    id:str = 'es_' + edit_date.replace('-','')
    ede_id:str = 'ese_' + edit_date.replace('-','')
    s = f'''{id} : createEditableStatement(deleted: false, editDate: \\"{edit_date}\\", field: \\"{field}\\", id: \\"{id}\\",statement: \\"{statement}\\"),'''
    s += f'{ede_id}: addEditableStatementEditor(editor:[\\"{editor_id}\\"], id:\\"{id}\\" ),'
    s += write_references(id,statement,pmid_extractor,reference_dict,journal_dict,author_dict)
    return s, id


def create_jax_description(id,field,description,editor_id,pmid_extractor:callable,reference_dict:dict,journal_dict:dict,author_dict:dict):
    s, es_id = createEditableStatement(description,field,editor_id,pmid_extractor,reference_dict,journal_dict,author_dict)
    s += f'addJaxGeneDescription(description:[\\"{es_id}\\"], id:\\"{id}\\"),'
    return s

def get_gene_id_from_entrez_id(entrez_id:str)->str:
    return 'geneInfo_gene_' + entrez_id


def get_omnigene_id_from_entrez_id(entrez_id:str)->str:
    return 'omnigene_' + entrez_id

def get_acessed_date_as_string(d:datetime)-> str:
    date_time = d.strftime("%m/%d/%Y")
    return date_time

def get_name_for_internet_reference(url:str ,accessed_date: str):
    pos1 = url.find('//') + 2
    pos2 = url.find('/',pos1)
    name = url[pos1:pos2]
    name += ' (accessed on:' + accessed_date + ')'
    return name

def create_myGeneInfo_gene(omni_gene:dict,auto_user_id,pmid_extractor:callable, reference_dict:dict, journal_dict:dict, author_dict:dict)->str:
    id = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    gene: str = omni_gene['symbol']
    chrom: str = omni_gene['chrom']
    strand: str = omni_gene['strand']
    start: int = omni_gene['start']
    end: int = omni_gene['end']
    entrez_id = omni_gene['entrez_gene_id']
    statement: str = gene + ' on chromosome ' + chrom + ' at ' + str(start) + '-' + str(end)
    if 'summary' in omni_gene:
        statement = omni_gene['summary']
    synonyms: str = '['
    for syn in omni_gene['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'
    s = f'{id}: createMyGeneInfo_Gene(chromosome: \\"{chrom}\\", end: {end}, entrezId: \\"{entrez_id}\\", id: \\"{id}\\", name: \\"{gene}\\" start: {start},  strand:{strand},synonyms: {synonyms}),'
    # statement:\\"{statement}\\",
    # addMyGeneInfo_GeneDescription(description: [ID!]!id: ID!): String
    field: str = 'geneDescription_' + id
    m, es_id= createEditableStatement(statement,field,auto_user_id,pmid_extractor,reference_dict,journal_dict,author_dict)
    s += m
    # addJaxGeneDescription(description: [ID!]!id: ID!): String
    s += f'addMyGeneInfo_GeneDescription(description:[\\"{es_id}\\"], id:\\"{id}\\"),'

    ref = omni_gene['reference']
    if ref['type'] == 'InternetReference':
        ref_id: str = 'ref_' + id
        accessed: str = get_acessed_date_as_string(ref['accessed_date'])
        ir_name: str = get_name_for_internet_reference(ref["url"], accessed)
        # createInternetReference(
        # accessed_date: String!
        # id: ID!
        # shortReference: String!
        # web_address: String!): String
        # Creates a InternetReference entity
        s += f'{ref_id}: createInternetReference(accessedDate:\\"{accessed}\\", id:\\"{ref_id}\\", shortReference: \\"{ir_name}\\", webAddress:\\"{ref["url"]}\\" ),'
        ref_id2 = 'gref_' + id
        s += f'{ref_id2}: addEditableStatementReferences(id:\\"{es_id}\\", references:[\\"{ref_id}\\"] ),'
    # print(s)
    return s

def create_uniprot_entry(omni_gene: dict, editor_id: str,pmid_extractor:callable,reference_dict:dict,journal_dict:dict,author_dict:dict)->str:
    mutation_payload: str = ''
    if 'sp_info' in omni_gene:
        sp_info = omni_gene['sp_info']
        id: str = sp_info['id']
        accessionNumber: str = sp_info['acc_num']
        statement: str = replace_characters(sp_info['function'])
        name: str = sp_info['name']
        uniprot_id: str = sp_info['uniprot_id']
        s = f'{id}: createUniprotEntry(accessionNumber: \\"{accessionNumber}\\", id: \\"{id}\\", name: \\"{name}\\",  uniprotID:\\"{uniprot_id}\\"),'
        mutation_payload += s
        # addUniprot_EntryGene(gene: [ID!]!id: ID!): String
        # Adds Gene to Uniprot_Entry entity
        gene_id = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
        s = f'addUniprotEntryGene(gene:[\\"{gene_id}\\"], id:\\"{id}\\" ),'
        mutation_payload += s
        field: str = 'proteinFunction_' + id
        m, es_id = createEditableStatement(statement, field, editor_id, pmid_extractor, reference_dict, journal_dict, author_dict)
        s += m
        # addUniprot_EntryFunction(function: [ID!]!id: ID!): String
        s += f'addUniprotEntryFunction(function:[\\"{es_id}\\"], id:\\"{id}\\" ),'
        # s += write_uniprot_references(es_id, statement, server)
        # s += write_references(es_id, statement, pmid_extractor, reference_dict, journal_dict, author_dict)
        mutation_payload += s
    return mutation_payload

def create_omniGene(omni_gene:dict, jax_gene_dict:dict, gene_description:str, editor_id, pmid_extractor:callable, reference_dict:dict, journal_dict:dict, author_dict:dict)->(str,str,str,str):
    id = get_omnigene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    gene: str = omni_gene['symbol']
    panel_name = omni_gene['panel_name']
    s = f'{id}: createOmniGene(id: \\"{id}\\", name: \\"{gene}\\", panelName:\\"{panel_name}\\" ),'
    # create geneDescription EditableStatement
    field1: str = 'geneDescription_' + id
    if gene_description==None:
        gene_description = '(Insert Gene Description)'
    statement1: str = gene_description
    (m, id1) = createEditableStatement(statement1,field1,editor_id,pmid_extractor, reference_dict, journal_dict, author_dict)
    s += m
#     addOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
#     Adds GeneDescription to OmniGene entity
    s += f'addOmniGeneGeneDescription(geneDescription:[\\"{id1}\\"], id:\\"{id}\\" ),'

    if 'category' in omni_gene and len(omni_gene['category'])>0:
        statement2: str = omni_gene['category']
    else:
        statement2: str = 'Neither'
        # create OncogenicCategory EditableStatement
    field2: str = 'OncogenicCategory_' + id
    (m, id2) = createEditableStatement(statement2,field2,editor_id,pmid_extractor, reference_dict, journal_dict, author_dict)
    s += m
     # addOmniGeneOncogenicCategory(id: ID!oncogenicCategory: [ID!]!): String
    # Adds OncogenicCategory to OmniGene entity
    s += f'addOmniGeneOncogenicCategory(id:\\"{id}\\", oncogenicCategory:[\\"{id2}\\"] ),'

    # addOmniGeneSynonymsString(id: ID!synonymsString: [ID!]!): String
    # Adds SynonymsString to OmniGene entity
    field3: str = 'SynonymsString_' + id
    statement3: str = gene
    if 'synonym' in omni_gene:
        statement3 = omni_gene['synonym']
    (m, id3) = createEditableStatement(statement3, field3, editor_id,pmid_extractor, reference_dict, journal_dict, author_dict)
    s += m
    s += f'addOmniGeneSynonymsString(id:\\"{id}\\", synonymsString:[\\"{id3}\\"] ),'

    # addOmniGeneJaxGene(id: ID!jaxGene: [ID!]!): String
# Adds JaxGene to OmniGene entity
#     jaxGene = get_gene_id_from_jax_id(omni_gene['entrez_gene_id'])
    if gene in jax_gene_dict:
        jaxGene = jax_gene_dict[gene]
        s += f'addOmniGeneJaxGene(id:\\"{id}\\", jaxGene:[\\"{jaxGene}\\"] ),'
    else:
        print("no jax gene for ",gene)
# addOmniGeneMyGeneInfoGene(id: ID!myGeneInfoGene: [ID!]!): String
# Adds MyGeneInfoGene to OmniGene entity
    myGeneInfoGene = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    s += f'addOmniGeneMyGeneInfoGene(id:\\"{id}\\", myGeneInfoGene:[\\"{myGeneInfoGene}\\"] ),'

    # addOmniGeneUniprot_entry(id: ID!uniprot_entry: [ID!]!): String
    # Adds Uniprot_entry to OmniGene entity
    if 'sp_info' in omni_gene:
        uniprot_id:str = omni_gene['sp_info']['id']
        s += f'addOmniGeneUniprotEntry(id:\\"{id}\\", uniprotEntry:[\\"{uniprot_id}\\"] ),'
    else:
        print("no Uniprot_entry for ", gene)

    return s, id, id2, id3


def return_graphql_boolean(b):
    s = 'false'
    if b:
        s = 'true'
    return s


# createUser(
# id: ID!
# isAdmin: Boolean!
# name: String!
# password: String!): String
# Creates a User entity
def write_users(users_dict:dict, server:str)->None:
    user_ids: dict = {}
    mutation_payload: str = ''
    for name, password in users_dict.items():
        id: str = 'user_' + datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
        user_ids[name] = id
        is_admin: bool = True
        mutation_payload += f'{id}: createUser(id:\\"{id}\\", isAdmin:{return_graphql_boolean(is_admin)}, name:\\"{name}\\", password:\\"{password}\\"),'
    # print(mutation_payload)
    send_mutation(mutation_payload, server)




def erase_neo4j(schema__graphql,server):
    uri = "bolt://" + server + ":7687"
    with open(schema__graphql, 'r') as file:
        idl_as_string = file.read()
    driver = GraphDatabase.driver(uri, auth=("neo4j", "omni"))
    with driver.session() as session:
        tx = session.begin_transaction()
        tx.run("match(a) detach delete(a)")
        result = tx.run("call graphql.idl('" + idl_as_string + "')")
        print(result.single()[0])
        tx.commit()
    driver.close()
