
import os
import re

FUSION_GENE_CAPTURE = "^\"[0-9]+?\"\s+\"" \
                      "(?P<gene1>.+?)" \
                      "\"\s+\".+?\"\s+\".+?\"\s+\".+?\"\s+\"" \
                      "(?P<gene2>.+?)\""

GENE_PDB_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+)," \
                   "(?P<pdbset1>\[?[',\s0-9A-Z]+\]?)," \
                   "(?P<gene2>[A-Z]+[0-9A-Z]+)," \
                   "(?P<pdbset2>\[?[',\s0-9A-Z]+\]?)"

FIINT_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+)," \
                "(?P<swissprot1>[A-Z]+[0-9A-Z]+)," \
                "(?P<gene2>[A-Z]+[0-9A-Z]+)," \
                "(?P<swissprot2>[A-Z]+[0-9A-Z]+)," \
                "(?P<interactome3d>.*)"

PROJ_DIR = "/Users/joshuaburkhart/Research/ReactomePPI"
OUT_DIR = PROJ_DIR + "/data/output"
FIINT_FN = OUT_DIR + "/FIInteract.txt"
ANNOT_PDB_FN = OUT_DIR + "/FIAnnotPdb.txt"
FUSION_RESULTS_FN = PROJ_DIR + "/gene-fusion-analysis/results/newDescription.txt"
FIINT_FUSION_INTERSECT_FN = OUT_DIR + "/FIInteractFusionEvents.txt"

fusion_gene_dict = dict()
gene_pdb_dict = dict()
fi_interactome_dict= dict()

def string2list(s):
    return s.replace('[', '') \
        .replace('\'', '') \
        .replace(',', '') \
        .replace(']', '') \
        .split()

# read FUSION_RESULTS_FN (if it exists)
if os.path.isfile(FUSION_RESULTS_FN):
    print("read FUSION_RESULTS_FN...")
    in_fptr = open(FUSION_RESULTS_FN)
    while 1:
        line = in_fptr.readline()
        if not line:
            break
        match = re.match(FUSION_GENE_CAPTURE, line)
        if match:
            # store map each gene to fusion event
            fusion_gene_dict[match.group('gene1')] = "{0}-{1}".format(match.group('gene1'),match.group('gene2'))
            fusion_gene_dict[match.group('gene2')] = "{0}-{1}".format(match.group('gene1'),match.group('gene2'))
    # close file
    in_fptr.close()

# open ANNOT_PDB_FN file
in_fptr = open(ANNOT_PDB_FN)
print("read ANNOT_PDB_FN...")
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(GENE_PDB_CAPTURE, line)
    if match:
        # store gene to pdb list map
        gene_pdb_dict[match.group('gene1')] = string2list(match.group('pdbset1'))
        gene_pdb_dict[match.group('gene2')] = string2list(match.group('pdbset2'))
# close file
in_fptr.close()

# open FIINT_FN file
in_fptr = open(FIINT_FN)
print("read FIINT_FN...")
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(FIINT_CAPTURE, line)
    if match:
        fi_interactome_dict[match.group('gene1')] = (match.group('gene1'),
                                   match.group('swissprot1'),
                                   match.group('gene2'),
                                   match.group('swissprot2'),
                                   match.group('interactome3d'))
        fi_interactome_dict[match.group('gene2')] = (match.group('gene1'),
                                   match.group('swissprot1'),
                                   match.group('gene2'),
                                   match.group('swissprot2'),
                                   match.group('interactome3d'))
# close file
in_fptr.close()

intersection = fusion_gene_dict.keys() & gene_pdb_dict.keys() & fi_interactome_dict.keys()
print("intersection yields {0} genes...".format(len(intersection)))

if len(intersection) > 0:
    with open(FIINT_FUSION_INTERSECT_FN, 'w') as out_fptr:
        out_fptr.write("FUSION EVENT,INTERACTION AFFECTED,"
                       "GENE 1 NAME, GENE 1 SWISSPROT, GENE 1 PDB SET,"
                       "GENE 2 NAME, GENE 2 SWISSPROT, GENE 2 PDB SET,"
                       "Interactome3D Result\n")
    # the thought is to later report only maximum pdb-pdb zdock scores (so as to fit them in one row) and link them to their zdock.out
    for gene in intersection:
        print("producing row for {0} with reported interaction {1}...".format(gene,fi_interactome_dict[gene]))
        with open(FIINT_FUSION_INTERSECT_FN, 'a') as out_fptr:
            out_fptr.write("F {0},I {1},{2},{3},{4},{5},{6},{7},{8}\n".format(
                fusion_gene_dict[gene].replace(',','~'),                              #FUSION EVENT
                "{0}-{1}".format(gene,fi_interactome_dict[gene][2]).replace(',','~'), #INTERACTION AFFECTED
                gene.replace(',','~'),                                                #GENE 1 NAME
                fi_interactome_dict[gene][1].replace(',','~'),                        #GENE 1 SWISSPROT
                '~'.join(map(str,gene_pdb_dict[gene])),                               #GENE 1 PDB SET
                fi_interactome_dict[gene][2].replace(',','~'),                        #GENE 2 NAME
                fi_interactome_dict[gene][3].replace(',','~'),                        #GENE 2 SWISSPROT
                '~'.join(map(str,gene_pdb_dict[fi_interactome_dict[gene][2]])),       #GENE 2 PDB SET
                fi_interactome_dict[gene][4].replace(',','~')                         #Interactome3D Result
            ))




