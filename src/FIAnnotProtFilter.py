import re
import mygene

GENE_DATA_ROW = '^(?P<gene1>[A-Z]+[0-9A-Z]+)\s+(?P<gene2>[A-Z]+[0-9A-Z]+)'

Gene_set = set()
FI_set = set()
Pdb_dict = dict()
AnnotProt_list = list()

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

[AnnotProt_list.append((gene_tup[0],
                        Pdb_dict[gene_tup[0]],
                        gene_tup[1],
                        Pdb_dict[gene_tup[1]]))
 for gene_tup in FI_set
 if gene_tup[0] in Pdb_dict and
    gene_tup[1] in Pdb_dict]

out_fptr = open("../data/output/FIAnnotProt.txt",'w+')
[out_fptr.write("{0},{1},{2},{3}\n".format(item[0],
                                           item[1],
                                           item[2],
                                           item[3]))
 for item in sorted(AnnotProt_list)]
out_fptr.close()
