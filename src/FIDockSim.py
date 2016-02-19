import os.path
import re
import itertools
from glob import glob
from ftplib import FTP
import gzip

PDB_CAPTURE = "^.+?," \
              "(?P<pdb1>[A-Z]+[0-9A-Z]+)," \
              ".+?," \
              "(?P<pdb2>[A-Z]+[0-9A-Z]+)"

GENE_PDB_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+)," \
                   "(?P<pdbset1>\[?[',\s0-9A-Z]+\]?)," \
                   "(?P<gene2>[A-Z]+[0-9A-Z]+)," \
                   "(?P<pdbset2>\[?[',\s0-9A-Z]+\]?)"

FIINT_CAPTURE = "^(?P<gene1>[A-Z]+[0-9A-Z]+)," \
                "(?P<swissprot1>[A-Z]+[0-9A-Z]+)," \
                "(?P<gene2>[A-Z]+[0-9A-Z]+)," \
                "(?P<swissprot2>[A-Z]+[0-9A-Z]+)," \
                "[^Interaction Reported by Interactome3D]"

ZDOCK_SCORE_CAPTURE = "^(-?[0-9]+.)?[0-9]*\t" \
                      "(-?[0-9]+.)?[0-9]*\t" \
                      "(-?[0-9]+.)?[0-9]*\t" \
                      "(-?[0-9]+.)?[0-9]*\t" \
                      "(-?[0-9]+.)?[0-9]*\t" \
                      "(?P<zdock_score>.*)"

prev_comp_pdb_set = set()
gene_pdb_dictionary = dict()
fi_no_interactome_set = set()


def string2list(s):
    return s.replace('[', '') \
        .replace('\'', '') \
        .replace(',', '') \
        .replace(']', '') \
        .split()


def download_pdb(pdb_id, download_dir):
    out_file = '{0}/{1}.pdb.gz'.format(download_dir, pdb_id)

    # see http://www.rcsb.org/pdb/static.do?p=download/ftp/index.html
    pdb_ftp = FTP("ftp.wwpdb.org")
    pdb_ftp.login()
    pdb_ftp.cwd("pub/pdb/data/structures/divided/pdb/{0}".format(pdb_id[1:3].tolower()))
    pdb_ftp.retrbinary('RETR pdb{0}.ent.gz'.format(pdb_id.tolower()),
                       open(out_file, 'wb').write)
    pdb_ftp.quit()
    os.system("gunzip {0}".format)


def zdock(pdb_receptor, pdb_ligand):
    os.system("mark_sur ../data/input/{0}.pdb ../data/input/{0}_m.pdb".format(pdb_receptor))
    os.system("mark_sur ../data/input/{0}.pdb ../data/input/{0}_m.pdb".format(pdb_ligand))
    os.system(
        "zdock -R ../data/input/{0}_m.pdb -L ../data/input/{1}_m.pdb -o ../data/input/zdock.out".format(pdb_receptor,
                                                                                                        pdb_ligand))


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

# write data/output/FIZdockScores.txt header: genes, pdbs, top zscore, number of complexes with top score
with open("../data/ouput/FIZdockScores.txt", 'w') as out_fptr:
    out_fptr.write("gene1,pdb_receptor,gene2,pdb_ligand,zdock_score,num_complexes\n")

# for stored FIInteract.txt rows
for interaction in fi_no_interactome_set:
    # for each pair of pdbs (some genes have more than 1) not previously computed
    gene1 = interaction[0]
    swissprot1 = interaction[1]
    gene2 = interaction[2]
    swissprot2 = interaction[3]
    for pdb_tup in itertools.product(gene_pdb_dictionary[gene1], gene_pdb_dictionary[gene2]):
        # check if .pdb file for each protein exists
        pdb1isfile = os.path.isfile("../data/input/{0}.pdb".format(pdb_tup[0]))
        pdb2isfile = os.path.isfile("../data/input/{0}.pdb".format(pdb_tup[1]))
        if (not pdb1isfile) and (not pdb2isfile):
            # delete any existing .pdb files & download .pdb files for each protein
            [os.remove(pdbfile) for pdbfile in glob('../data/input/*.pdb')]
            download_pdb(pdb_tup[0], "../data/input")
            download_pdb(pdb_tup[1], "../data/input")
        elif pdb1isfile:
            # download .pdb 2
            download_pdb(pdb_tup[1], "../data/input")
        elif pdb2isfile:
            # download .pdb 1
            download_pdb(pdb_tup[0], "../data/input")

        # perform zdock execution on .pdb files
        zdock(pdb_tup[0], pdb_tup[1])

        # record top zdock score and number of rows (complexes) with that score
        top_zdock_score = float("-inf")
        num_top_scores = 0
        in_fptr = open("../data/input/zdock.out")
        while 1:
            line = in_fptr.readline()
            if not line:
                break
            match = re.match(ZDOCK_SCORE_CAPTURE, line)
            if match and float(match.group('zdock_score')) >= top_zdock_score:
                top_zdock_score = float(match.group('zdock_score'))
                num_top_scores += 1
            else:
                break  # zdock scores are listed hi -> lo
        in_fptr.close()

        # open data/output/FIZdockScores.txt
        # store genes, pdbs, top zscore, number of complexes with top score
        with open("../data/ouput/FIZdockScores.txt", 'a') as out_fptr:
            out_fptr.write("{0},{1},{2},{3},{4},{5}\n".format(gene1, pdb_tup[0],
                                                              gene2, pdb_tup[1],
                                                              top_zdock_score, num_top_scores))
