import re
import mygene
import urllib.request
import urllib.parse
import xml

GENE_DATA_ROW = '^(?P<gene1>[A-Z]+[0-9A-Z]+)\s+(?P<gene2>[A-Z]+[0-9A-Z]+)'

Gene_set = set()
FI_set = set()
Pdb_dict = dict()
Uni_dict = dict()
AnnotPdb_list = list()
AnnotUni_list = list()

in_fptr = open("../data/input/FIsInGene_121514_with_annotations.txt")
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(GENE_DATA_ROW, line)
    if match:
        gene_tup = (match.group('gene1'), match.group('gene2'))
        Gene_set.add(gene_tup[0])
        Gene_set.add(gene_tup[1])
        FI_set.add(gene_tup)
in_fptr.close()

mg = mygene.MyGeneInfo()
annot_gene_list = mg.querymany(Gene_set,
                               scopes="symbol",
                               fields=["pdb", "uniprot"],
                               species="human",
                               as_dataframe=False)

for gene in annot_gene_list:
    try:
        Pdb_dict[gene['query']] = gene['pdb']  # this can return more than one pdb
    except KeyError:
        # ignore genes without PDB's
        pass
    try:
        Uni_dict[gene['query']] = gene['uniprot']['Swiss-Prot']
    except KeyError:
        # ignore genes without uniprot's
        pass


[AnnotPdb_list.append((gene_tup[0],
                        Pdb_dict[gene_tup[0]],
                        gene_tup[1],
                        Pdb_dict[gene_tup[1]]))
 for gene_tup in FI_set
 if gene_tup[0] in Pdb_dict and
    gene_tup[1] in Pdb_dict]

[AnnotUni_list.append((gene_tup[0],
                        Uni_dict[gene_tup[0]],
                        gene_tup[1],
                        Uni_dict[gene_tup[1]]))
 for gene_tup in FI_set
 if gene_tup[0] in Uni_dict and
    gene_tup[1] in Uni_dict]

out_fptr = open("../data/output/FIAnnotPdb.txt",'w+')
[out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                           item[1],
                                           item[2],
                                           item[3]))
 for item in sorted(AnnotPdb_list)]
out_fptr.close()

out_fptr = open("../data/output/FIAnnotUni.txt",'w+')
[out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                           item[1],
                                           item[2],
                                           item[3]))
 for item in sorted(AnnotUni_list)]
out_fptr.close()

#not yet implemented, bail out
exit(1)
for line in AnnotUni_list:
    gene1 = line[0]
    unip1 = line[1]
    gene2 = line[2]
    unip2 = line[3]

    query_data = {}
    query_data['queryProt1'] = unip1
    query_data['queryProt2'] = unip2

    url_values = urllib.parse.urlencode(query_data)

    url = 'http://interactome3d.irbbarcelona.org/api/getInteractionStructures'
    full_url = '{0}?{1}'.format(url,url_values)

    #TODO: Parse the XML
    #https://docs.python.org/3/library/xml.etree.elementtree.html
    #http://interactome3d.irbbarcelona.org/help.php#restful
    #http://interactome3d.irbbarcelona.org/help.php#interactions_dat_file

    with urllib.request.urlopen(full_url) as response:
        html = response.read()
        xmlTree = xml.etree.ElementTree.parse(html)
        xmlRoot = root = xmlTree.getroot()
