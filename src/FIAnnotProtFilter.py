import re
import mygene

GENE_DATA_ROW = '^(?P<gene1>[A-Z]+[0-9A-Z]+)\s+(?P<gene2>[A-Z]+[0-9A-Z]+)'

Gene_set = set.new
FI_set = set.new
Pdb_dict = hash.new
AnnotProt_list = list.new

in_fptr = open("../data/input/FIsInGene_121514_with_annotations.txt")
while 1:
    line = in_fptr.readline()
    if not line:
        break;
    if re.match(GENE_DATA_ROW, line):
        gene_tup = (re.group('gene1'),re.group('gene2'))
        Gene_set.add(gene_tup[0])
        Gene_set.add(gene_tup[1])
        FI_set.add(gene_tup)
in_fptr.close()

mg = mygene.MyGeneInfo()
annot_gene_list = mg.querymany(Gene_set,
             scopes="symbol",
             fields=["pdb","uniprot"],
             species="human",
             as_dataframe=False)

for gene in annot_gene_list:
    try:
        Pdb_dict[gene['query']] = gene['pdb'] # this can return more than one pdb
    except KeyError:
        # ignore genes without PDB's
        pass

for gene_tup in FI_set:
    if gene_tup[0] in Pdb_dict and gene_tup[1] in Pdb_dict:
        AnnotProt_list.add((gene_tup[0],Pdb_dict[0],gene_tup[1],Pdb_dict[1]))

out_fptr = open("../data/output/FIAnnotProt.txt")
for item in AnnotProt_list:
    out_fptr.write(item[0],item[1],item[2],item[3])
out_fptr.close()
