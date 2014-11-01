__author__ = 'pshort'

from Bio import Entrez
from lxml import etree

#Using biopython to get pubmed metadata from NCBI

def get_pubmed_metadata(pmid_list, email):
    """
    use biopython.entrez to return pubmed metadata as a dict
    """

    post_size = len(pmid_list)
    print('Preparing to post %i requests to PubMed.' % post_size)
    Entrez.email = email
    Entrez.tool = 'biopython'
    search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(pmid_list)))
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    batch_size = 500
    total_words = 0
    punc = [",", ".", "/", ";", "'", "?", "&", "-", ")", "(", "\"", "\\"]

    pmid2abstract_data = {}
    pmid2authors_data = {}
    pmid2affil_data = {}
    pmid2mesh_data = {}
    pmid2journal_data = {}
    pmid2print_medium = {}
    for start in range(0, post_size, batch_size):
        out_handle = open('pubmed_abstract_test.xml', 'w')
        end = min(post_size, start + batch_size)
        print("Going to download record %i to %i" % (start + 1, end))
        fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="xml", retstart=start, retmax=batch_size,
                                     webenv=webenv, query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
        out_handle.close()
        xml_data = 'pubmed_abstract_test.xml'
        tree = etree.parse(xml_data)
        root = tree.getroot()

        for article in root:
            pmid = article.find('.//PMID').text
            abstract = {}
            try:
                abtext = article.find('.//AbstractText').text
                abtext = abtext.split(" ")
                total_words += len(abtext)
                for word in abtext:
                    for p in punc:
                        word = word.replace(p, "")
                    if word not in abstract:
                        abstract[word] = 1
                pmid2abstract_data[pmid] = abstract.keys()
            except AttributeError:
                pmid2abstract_data[pmid] = ['None']
                pass
            authors = {}
            try:
                author_total = article.find('.//AuthorList')
                for author in author_total:
                    lastname = author[0].text
                    if len(author) > 1:
                        firstname = author[1].text
                        name = unicode(firstname) + " " + unicode(lastname)
                    else:
                        name = lastname
                    name = name.replace("\"", "")
                    if name not in authors:
                        authors[name] = 1
                pmid2authors_data[pmid] = authors.keys()
            except:
                pmid2authors_data[pmid] = ['None']
                pass
            try:
                affiliation = article.find('.//Affiliation').text
                affiliation_list = affiliation.split(",")
                if len(affiliation_list) > 1:
                    affiliation = affiliation_list[0] + ", " + affiliation_list[1]
                pmid2affil_data[pmid] = affiliation
            except AttributeError:
                pmid2affil_data[pmid] = 'None'
                pass
            mesh_headings = {}
            try:
                mesh_list = article.find('.//MeshHeadingList')
                for mesh in mesh_list:
                    mesh_descriptor = mesh[0].text
                    if mesh_descriptor not in mesh_headings:
                        mesh_headings[mesh_descriptor] = 1
                pmid2mesh_data[pmid] = mesh_headings.keys()
            except AttributeError and TypeError:
                pmid2mesh_data[pmid] = ['None']
                pass
            try:
                title = article.find('.//Title').text
                pmid2journal_data[pmid] = title
            except AttributeError:
                pmid2journal_data[pmid] = 'None'
                pass
            try:
                medium = article.find('.//JournalIssue').get('CitedMedium')
                pmid2print_medium[pmid] = medium
            except AttributeError:
                pmid2print_medium[pmid] = 'None'
                pass

    pmid2metadata = {}
    for pmid in pmid2authors_data:
        data_dict = {'author': pmid2authors_data[pmid], 'journal': pmid2journal_data[pmid],
                     'affiliation': pmid2affil_data[pmid],
                     'abs_word': pmid2abstract_data[pmid], 'mesh_head': pmid2mesh_data[pmid],
                     'print_med': pmid2print_medium[pmid]}
        pmid2metadata[pmid] = data_dict

    return pmid2metadata