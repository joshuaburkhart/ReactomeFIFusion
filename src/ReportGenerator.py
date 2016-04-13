from lxml import html
import requests
import urllib
import os

# Parse Input

# Loop over Interactions

# Interactome3D

page = requests.get('http://interactome3d.irbbarcelona.org/interaction.php?ids=P63104;Q02156&dataset=human&rs=True&connect=1')
tree = html.fromstring(page.content)

#Interaction Complex Name
intcn_cplx_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[5]/text()'

#Complex 1: chain (B)
cplx1_chain_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[6]/text()'

#Complex 1: length (1 P63104 245 )
cplx1_start_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[1]/div/table/tr/td[1]/text()'
cplx1_end_xpath   = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[1]/div/table/tr/td[3]/text()'

#Complex 1: from/to coordinates for chain (2      228)
cplx1_from_xpath  = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[8]/text()'
cplx1_to_xpath    = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[9]/text()'

#Complex 2: chain (Q)
cplx2_chain_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[12]/text()'

#Complex 2: length (1 Q02156 737)
cplx2_start_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[2]/div/table/tr/td[1]/text()'
cplx2_end_xpath   = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[2]/div/table/tr/td[3]/text()'

#Complex 2: from/to coordinates for chain (365      372)
cplx2_from_xpath  = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[14]/text()'
cplx2_to_xpath    = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[15]/text()'

intcn_cplx_val  = tree.xpath(intcn_cplx_xpath)[0]
cplx1_chain_val = tree.xpath(cplx1_chain_xpath)[0]
cplx1_start_val = tree.xpath(cplx1_start_xpath)[0]
cplx1_end_val   = tree.xpath(cplx1_end_xpath)[0]
cplx1_from_val  = tree.xpath(cplx1_from_xpath)[0]
cplx1_to_val    = tree.xpath(cplx1_to_xpath)[0]
cplx2_chain_val = tree.xpath(cplx2_chain_xpath)[0]
cplx2_start_val = tree.xpath(cplx2_start_xpath)[0]
cplx2_end_val   = tree.xpath(cplx2_end_xpath)[0]
cplx2_from_val  = tree.xpath(cplx2_from_xpath)[0]
cplx2_to_val    = tree.xpath(cplx2_to_xpath)[0]

print('Interaction Complex: {0}'.format(intcn_cplx_val))
print('Complex 1 Chain: {0}'.format(cplx1_chain_val))
print('Complex 1 Length: {0} - {1}'.format(cplx1_start_val,cplx1_end_val))
print('Complex 1 Included in Interaction: {0} - {1}'.format(cplx1_from_val,cplx1_to_val))
print('Complex 2 Chain: {0}'.format(cplx2_chain_val))
print('Complex 2 Length: {0} - {1}'.format(cplx2_start_val,cplx2_end_val))
print('Complex 2 Included in Interaction: {0} - {1}'.format(cplx2_from_val,cplx2_to_val))

# Create Directories

PDB_DIR = '{0}/../data/output/reports/{1}/src/pdb'.format(os.path.dirname(os.path.realpath(__file__)),intcn_cplx_val)
PNG_DIR = '{0}/../data/output/reports/{1}/src/png'.format(os.path.dirname(os.path.realpath(__file__)),intcn_cplx_val)
MD_DIR = '{0}/../data/output/reports/{1}/md'.format(os.path.dirname(os.path.realpath(__file__)),intcn_cplx_val)
HTML_DIR = '{0}/../data/output/reports/{1}/html'.format(os.path.dirname(os.path.realpath(__file__)),intcn_cplx_val)
PDF_DIR = '{0}/../data/output/reports/{1}/pdf'.format(os.path.dirname(os.path.realpath(__file__)),intcn_cplx_val)

if not os.path.exists(PDB_DIR):
    os.makedirs(PDB_DIR)
if not os.path.exists(PNG_DIR):
    os.makedirs(PNG_DIR)
if not os.path.exists(MD_DIR):
    os.makedirs(MD_DIR)
if not os.path.exists(HTML_DIR):
    os.makedirs(HTML_DIR)
if not os.path.exists(PDF_DIR):
    os.makedirs(PDF_DIR)

# Download PDB
pdb_gz_path = '{0}/{1}.pdb.gz'.format(PDB_DIR,intcn_cplx_val)
if os.path.exists(pdb_gz_path):
    os.remove(pdb_gz_path)
urllib.request.urlretrieve('http://files.rcsb.org/download/{0}.pdb.gz'.format(intcn_cplx_val.upper()), pdb_gz_path)
pdb_path = '{0}/{1}.pdb'.format(PDB_DIR,intcn_cplx_val)
if os.path.exists(pdb_path):
    os.remove(pdb_path)
os.system('gunzip {0}'.format(pdb_gz_path))

# Generate Image
png_path = '{0}/{1}.png'.format(PNG_DIR,intcn_cplx_val)
if os.path.exists(png_path):
    os.remove(png_path)
os.system('java -Djava.awt.headless=false -Xmx512m -jar "/Applications/Jmol.jar" -onj \
"load file {0}; \
select all; \
color grey; \
select:{1}; \
color yellow; \
select:{2}; \
color aqua; \
show BEST ROTATION; \
write {3}"'.format(pdb_path,cplx1_chain_val,cplx2_chain_val,png_path))

# PDBE

pdbe_root_url = 'http://www.ebi.ac.uk/pdbe/entry/pdb/{0}'.format(intcn_cplx_val)

#Biology
page = requests.get('{0}/{1}'.format(pdbe_root_url,'biology'))
tree = html.fromstring(page.content)

rctn_xpath = '/html/body/div[2]/div/section/div[3]/div/section[1]/div/div/section[1]/div/text()'

rctn_val = tree.xpath(rctn_xpath)[0].strip()

print('Reaction: {0}'.format(rctn_val))

#Experiment
page = requests.get('{0}/{1}'.format(pdbe_root_url,'experiment'))
tree = html.fromstring(page.content)

dscrp_xpath = '/html/body/div[2]/div/section/div[1]/div[3]/div/section[1]/div/div[1]/div/text()'

dscrp_val = tree.xpath(dscrp_xpath)[0].strip()

print('Author Description: {0}'.format(dscrp_val))

#Citation
page = requests.get('{0}/{1}'.format(pdbe_root_url,'citations'))
tree = html.fromstring(page.content)

cit_xpath = '/html/body/div[2]/div/section/div[1]/div[2]/div[2]/div/div[3]/div[2]/a/@href'

cit_val = tree.xpath(cit_xpath)[0].strip()

print('Citation: {0}'.format(cit_val))