# coding=utf-8
import re
import mygene
import urllib.request
import urllib.parse
import xml.etree.ElementTree

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

# Print Genes With PDB's to File
out_fptr = open("../data/output/FIAnnotPdb.txt", 'w+')
[out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                           item[1],
                                           item[2],
                                           item[3]))
 for item in sorted(AnnotPdb_list)]
out_fptr.close()

# Print Genes With UNIPROT: Swiss-Prot's to File
out_fptr = open("../data/output/FIAnnotUni.txt", 'w+')
[out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                           item[1],
                                           item[2],
                                           item[3]))
 for item in sorted(AnnotUni_list)]
out_fptr.close()

# Print Genes With PDB's to File With Their UNIPROT: Swiss-Prot's
out_fptr = open("../data/output/FIAnnotPdbWithUni.txt", 'w+')
for item in sorted(AnnotPdb_list):
    try:
        out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                                  Uni_dict[item[0]],
                                                  item[2],
                                                  Uni_dict[item[2]]))
    except KeyError:
        # ignore genes without both PDB's and Swiss-Prot's
        pass
out_fptr.close()

intr_count = 0
non_intr_count = 0
except_count = 0
loop_count = 0
loop_len = len(AnnotPdb_list)
out_fptr = open("../data/output/FIInteract.txt", 'w+')
for line in sorted(AnnotPdb_list):
    try:
        loop_count += 1
        gene1 = line[0]
        unip1 = Uni_dict[line[0]]
        gene2 = line[2]
        unip2 = Uni_dict[line[2]]

        query_data = {}
        query_data['queryProt1'] = unip1
        query_data['queryProt2'] = unip2

        url_values = urllib.parse.urlencode(query_data)

        url = 'http://interactome3d.irbbarcelona.org/api/getInteractionStructures'
        full_url = '{0}?{1}'.format(url, url_values)

        print("Querying interaction between {0} and {1} ({2} of {3})...".format(
            gene1,
            gene2,
            loop_count,
            loop_len
        ))

        try:
            with urllib.request.urlopen(full_url) as response:
                html = response.read()

            xmlTree = xml.etree.ElementTree.fromstring(html)
            if len(xmlTree.getchildren()[0].getchildren()) > 0:
                intr_msg = 'Interaction Reported by Interactome3D,{0}'.format(
                    full_url
                )
                intr_count += 1
            else:
                intr_msg = 'No Interaction Reported by Interactome3D,'
                non_intr_count += 1
        except Exception as e:
            intr_msg = 'Exception Raised: {0},'.format(e)
            except_count += 1
        out_fptr.write("{0},{1},{2},{3},{4}\n".format(
            gene1,
            unip1,
            gene2,
            unip2,
            intr_msg
        ))
    except KeyError:
        # Ignore genes w/o PDB's & Swiss-Prots
        pass
out_fptr.close()
print('Interaction Count: {0}\n' \
      'Non Interaction Count: {1}\n' \
      'Exception Count: {2}'.format(
    intr_count,
    non_intr_count,
    except_count
))
