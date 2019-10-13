from bs4 import BeautifulSoup
import glob
import re

def get_addgene_dict(f):
    with open(f, 'r') as fh:
        html_doc = fh.read()
    soup = BeautifulSoup(html_doc, 'html.parser') 

    name = soup.title.contents[0] # Get Name
    url = soup.head.find('link',{'rel':'canonical'})['href'] # Get url
    num = re.findall('.org\/(.*?)\/',url)
    try:
        purpose = soup.find('div', {"class": "field-content"}).contents[0].strip() # Get purpose
    except:
        purpose = None
    description = soup.head.find('meta', {'name': 'description'})['content'] # Get description

    dct = {'name':name,'num':num[0],'purpose':purpose,'description':description,'backbone_size':None,'Selectable markers':None,'Alt name':None,'Entrez Gene':None,'tags':[],'Insert Size (bp)':None,
            'Backbone manufacturer':None,'Species':None,'Promoter':None,'5′ sequencing primer':None, 'Cloning method':None,'Gene/Insert name':None,'Vector type':None, 'Vector backbone':None,
            'Bacterial Resistance(s)':None,'Growth Temperature':None,'Growth Strain(s)':None,'Copy number':None,'Terms and Licenses':[]}

    for x in soup.findAll('li', {'class': 'field'}):
        d = re.findall('<\/div>(.*?)<\/li>', str(x), re.DOTALL)
        if d == []:
            d = re.findall('<\/span>(.*?)<\/li>', str(x), re.DOTALL)
        if d != []:
            #print(d)
            dm = d[0].strip()
            div = x.findAll('div')
            if div == []:
                div = x.findAll('span')
            div = div[0].contents[0]
            if div == 'Vector backbone':
                vec_backbone = dm.split('\n')[0].strip()
                #print('Vector Backbone: {}'.format(vec_backbone))
                dct['Vector backbone'] = vec_backbone

            if div == 'Backbone manufacturer':
                backbone_man = dm
                #print('Backbone manufacturer: {}'.format(backbone_man))
                dct['Backbone manufacturer'] = backbone_man

            backbone_size_str = """
                                Backbone size

                                w/o insert

                                (bp)"""
            if div == backbone_size_str:
                backbone_size = int(dm)
                #print('Backbone size: {}'.format(backbone_size))
                dct['backbone_size'] = backbone_size

            if div == 'Vector type':
                vec_type = dm
                #print('Vector type: {}'.format(vec_type))
                dct['Vector type'] = vec_type

            if div == 'Selectable markers':
                select_marker = dm
                #print('Selectable markers: {}'.format(select_marker))
                dct['Selectable markers'] = select_marker

            if div == 'Bacterial Resistance(s)':
                b_res = dm
                #print('Bacterial Resistance(s): {}'.format(b_res))
                dct['Bacterial Resistance(s)'] = b_res

            if div == 'Growth Temperature':
                growth_temp = dm
                #print('Growth Temperature: {}'.format(growth_temp))
                dct['Growth Temperature'] = growth_temp

            if div == 'Growth Strain(s)':
                strain = dm
                #print('Growth Strain(s): {}'.format(strain))
                dct['Growth Strain(s)'] = strain

            if div == 'Copy number':
                copy_num = dm
                #print('Copy number: {}'.format(copy_num))
                dct['Copy number'] = copy_num

            if div == 'Gene/Insert name':
                gene_insert = dm
                #print('Gene/Insert name: {}'.format(gene_insert))
                dct['Gene/Insert name'] = gene_insert

            if div == 'Alt name':
                alt_name = dm
                #print('Alt name: {}'.format(alt_name))
                dct['Alt name'] = alt_name

            if div == 'Species':
                species = dm
                #print('Species: {}'.format(species))
                dct['Species'] = species

            if div == 'Insert Size (bp)':
                insert_size = int(dm)
                #print('Insert Size (bp): {}'.format(insert_size))
                dct['Insert Size (bp)'] = insert_size

            if div == 'Entrez Gene':
                a = re.findall('<a href="(.*?)"',dm)
                entrez = a[0]
                #print('Entrez Gene: {}'.format(entrez))
                dct['Entrez Gene'] = entrez

            if div == 'Terms and Licenses':
                d = re.findall('<\/div>(.*?)<\/ul>', str(x), re.DOTALL)[0]
                a = re.findall('">(.*?)<\/a><\/li>',d)
                terms_license = a
                #print('Terms and Licenses: {}'.format(terms_license))
                dct['Terms and Licenses'] = terms_license

            if div == 'Promoter':
                promoter = dm
                #print('Promoter: {}'.format(promoter))
                dct['Promoter'] = promoter

            tag_str = """
            Tag
            / Fusion Protein"""
            if div == tag_str:
                d = re.findall('<\/span>(.*?)<\/ul>', str(x), re.DOTALL)[0]
                a = re.findall('<li>(.*?)<\/li>',d)
                tags = a
                #print('Tags: {}'.format(tags))
                dct['tags'] = tags

            if div == 'Cloning method':
                cloning_method = dm
                #print('Cloning method: {}'.format(cloning_method))
                dct['Cloning method'] = cloning_method

            if div == '5′ sequencing primer':
                seq_primer = dm
                #print('Sequencing primer: {}'.format(seq_primer))
                dct['5′ sequencing primer'] = seq_primer
                
    return dct

import sqlite3

conn = sqlite3.connect('test.db')
cursor=conn.cursor()
table_def = """CREATE TABLE plasmids
(
id INTEGER PRIMARY KEY NOT NULL,
num text NOT NULL,
name text NOT NULL,
purpose text,
description text,
vector_backbone text,
backbone_manufacturer text,
backbone_size int,
vector_type text,
selectable_marker text,
bacterial_resistance text,
growth_temperature text,
growth_strain text,
copy_number text,
gene_name text,
alt_name text,
species text,
insert_size int,
entrez text,
promoter text,
cloning_method text,
seq_primer text
);

CREATE TABLE tags
(
id INTEGER PRIMARY KEY NOT NULL,
tag text UNIQUE
);

CREATE TABLE terms
(
id INTEGER PRIMARY KEY NOT NULL,
term text UNIQUE

);

CREATE TABLE plasmid_tags
(
id INTEGER PRIMARY KEY NOT NULL,
PlasmidId INTEGER NOT NULL,
TagId INTEGER NOT NULL,

FOREIGN KEY (PlasmidId) REFERENCES plasmids(id),
FOREIGN KEY (TagId) REFERENCES tags(id)
);

CREATE TABLE plasmid_terms
(
id INTEGER PRIMARY KEY NOT NULL,
PlasmidId INTEGER NOT NULL,
TermId INTEGER NOT NULL,

FOREIGN KEY (PlasmidId) REFERENCES plasmids(id),
FOREIGN KEY (TermId) REFERENCES terms(id)
)"""
for x in table_def.split(';'):
    print('=')
    conn.execute(x)

def insert_plasmid_string(dct,conn=conn):
    cursor=conn.cursor()
    plasmid_insert = """INSERT INTO  plasmids (num,name,purpose,description,vector_backbone,backbone_manufacturer,backbone_size,
    vector_type, selectable_marker, bacterial_resistance, growth_temperature, growth_strain, copy_number,
    gene_name, alt_name, species, insert_size, entrez, promoter, cloning_method, seq_primer) 
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
    
    
    
    cursor.execute(plasmid_insert,(dct['num'],dct['name'],dct['purpose'],dct['description'],dct['Vector backbone'],dct['Backbone manufacturer'], dct['backbone_size'],dct['Vector type'], dct['Selectable markers'],dct['Bacterial Resistance(s)'],dct['Growth Temperature'],dct['Growth Strain(s)'],dct['Copy number'], dct['Gene/Insert name'], dct['Alt name'], dct['Species'],dct['Insert Size (bp)'], dct['Entrez Gene'], dct['Promoter'], dct['Cloning method'], dct['5′ sequencing primer']))
    
    plasmid_id = cursor.execute("SELECT id FROM plasmids WHERE num = '{}'".format(dct['num'])).fetchone()[0]
    
    #term_insert = "INSERT OR IGNORE INTO terms (term) VALUES ({})"
    
    for tag in dct['tags']:
        cursor.execute("INSERT OR IGNORE INTO tags (tag) VALUES ('{}')".format(tag))
        tag_id = cursor.execute("SELECT id FROM tags WHERE tag = '{}'".format(tag)).fetchone()[0]
        cursor.execute("INSERT INTO plasmid_tags (PlasmidId,TagId) VALUES ('{}','{}')".format(plasmid_id,tag_id)) 
    
    for term in dct['Terms and Licenses']:
        cursor.execute("INSERT OR IGNORE INTO terms (term) VALUES ('{}')".format(term))
        term_id = cursor.execute("SELECT id FROM terms WHERE term = '{}'".format(term)).fetchone()[0]
        cursor.execute("INSERT INTO plasmid_terms (PlasmidId,TermId) VALUES ('{}','{}')".format(plasmid_id,term_id)) 
    
    conn.commit()
    return ('Upload of {}'.format(dct['num']))


for f in glob.glob('./../gene_list/*'):
    insert_plasmid_string(get_addgene_dict(f), conn=conn)
