import os.path
import re

PDB_CAPTURE = "^.+?,(?P<pdb1>[A-Z]+[0-9A-Z]+)," \
              ".+?,(?P<pdb2>[A-Z]+[0-9A-Z]+)"
GENE_PDB_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+),(?P<pdbset1>\[?[',\s0-9A-Z]+\]?)," \
                   "(?P<gene2>[A-Z]+[0-9A-Z]+),(?P<pdbset2>\[?[',\s0-9A-Z]+\]?)"
FIINT_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+),(?P<swissprot1>[A-Z]+[0-9A-Z]+)," \
                "(?P<gene2>[A-Z]+[0-9A-Z]+),(?P<swissprot2>[A-Z]+[0-9A-Z]+)," \
                "[^Interaction Reported by Interactome3D]"

prev_comp_pdb_set = set()
gene_pdb_dictionary = dict()
fi_no_interactome_set = set()


def string2list(s):
    return s.replace('[','') \
        .replace('\'','') \
        .replace(',','') \
        .replace(']','') \
        .split()

# read data/output/FIZdockScores (if it exists)
fi_zdock_scores_path = "../data/input/FIZdockScores.txt"
if os.path.isfile(fi_zdock_scores_path):
    in_fptr = open(fi_zdock_scores_path)
    while 1:
        line = in_fptr.readline()
        if not line:
            break
        match = re.match(PDB_CAPTURE, line)
        if match:
            # store previously computed pdb pairs
            pdb_tup = (match.group('pdb1'), match.group('pdb2'))
            prev_comp_pdb_set.add(pdb_tup)
    # close file
    in_fptr.close()

# open data/output/FIAnnotPdb.txt file
in_fptr = open("../data/input/FIAnnotPdb.txt")
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(GENE_PDB_CAPTURE, line)
    if match:
        # store gene to pdb list map
        gene_pdb_dictionary[match.group('gene1')] = string2list(match.group('pdbset1'))
        gene_pdb_dictionary[match.group('gene2')] = string2list(match.group('pdbset2'))
# close file
in_fptr.close()

# open data/output/FIInteract.txt file
in_fptr = open("../data/input/FIInteract.txt")
while 1:
    line = in_fptr.readline()
    if not line:
        break
    match = re.match(FIINT_CAPTURE, line)
    if match:
        # store rows without Interactome3D results in memory
        fi_no_interactome_set.add((match.group('gene1'),
                                   match.group('swissprot1'),
                                   match.group('gene2'),
                                   match.group('swissprot2')))
# close file
in_fptr.close()

# for stored FIInteract.txt rows

# for each pair of pdbs (some genes have more than 1) not previously computed

# check if .pdb file for each protein exists

# if not, delete any existing .pdb files & download .pdb files for each protein

# perform zdock execution on .pdb files

# record top zdock score and number of rows (complexes) with that score

# open data/output/FIZdockScores.txt for writing
# store genes, pdbs, top zscore, number of complexes with top score
# close file
