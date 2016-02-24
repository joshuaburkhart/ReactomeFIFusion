import ftplib
import os.path
import re
import itertools
from glob import glob
from ftplib import FTP

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

OUT_DIR = "/home/burkhart/Software/ReactomePPI/data/output"
FIINT_FN = OUT_DIR + "/FIInteract.txt"
ANNOT_PDB_FN = OUT_DIR + "/FIAnnotPdb.txt"
ZDOCK_OUT_FN = OUT_DIR + "/zdock.out"
ZDOCK_SCORE_FN = OUT_DIR + "/FIZdockScores.txt"

ZDOCK_DIR = "/home/burkhart/zdock3.0.2_linux_x64"
ZDOCK = ZDOCK_DIR + "/zdock"
MARK_SUR = ZDOCK_DIR + "/mark_sur"

prev_comp_pdb_set = set()
gene_pdb_dictionary = dict()
fi_no_interactome_set = set()


def string2list(s):
    print("string2list...")
    return s.replace('[', '') \
        .replace('\'', '') \
        .replace(',', '') \
        .replace(']', '') \
        .split()


def download_pdb(pdb_id, download_dir):
    print("download_pdb...")
    out_file = '{0}/{1}.pdb.gz'.format(download_dir, pdb_id)

    # see http://www.rcsb.org/pdb/static.do?p=download/ftp/index.html
    pdb_ftp = FTP("ftp.wwpdb.org")
    pdb_ftp.login()
    pdb_ftp.cwd("pub/pdb/data/structures/divided/pdb/{0}".format(pdb_id[1:3].lower()))
    pdb_ftp.retrbinary('RETR pdb{0}.ent.gz'.format(pdb_id.lower()),
                       open(out_file, 'wb').write)
    pdb_ftp.quit()
    os.system("gunzip -f {0}".format(out_file))


def zdock(pdb_receptor, pdb_ligand):
    print("zdock...")

    cp_R_to_zdock_cmd = "cp {0}/{1}.pdb {2}/".format(OUT_DIR, pdb_receptor, ZDOCK_DIR)
    cp_L_to_zdock_cmd = "cp {0}/{1}.pdb {2}/".format(OUT_DIR, pdb_ligand, ZDOCK_DIR)
    os.system(cp_R_to_zdock_cmd)
    os.system(cp_L_to_zdock_cmd)

    mark_sur_R_cmd = "cd {0} && {1} {2}.pdb {2}_m.pdb".format(ZDOCK_DIR, MARK_SUR, pdb_receptor)
    mark_sur_L_cmd = "cd {0} && {1} {2}.pdb {2}_m.pdb".format(ZDOCK_DIR, MARK_SUR, pdb_ligand)
    print(mark_sur_R_cmd)
    os.system(mark_sur_R_cmd)
    print(mark_sur_L_cmd)
    os.system(mark_sur_L_cmd)

    zdock_cmd = "cd {0} && {1} -R {2}_m.pdb -L {3}_m.pdb -o {4}".format(ZDOCK_DIR, ZDOCK, pdb_receptor, pdb_ligand,
                                                                        ZDOCK_OUT_FN)
    print(zdock_cmd)
    os.system(zdock_cmd)

    cp_R_to_data_cmd = "cp {0}/{1}_m.pdb {2}/".format(ZDOCK_DIR, pdb_receptor, OUT_DIR)
    cp_L_to_data_cmd = "cp {0}/{1}_m.pdb {2}/".format(ZDOCK_DIR, pdb_ligand, OUT_DIR)
    os.system(cp_R_to_data_cmd)
    os.system(cp_L_to_data_cmd)

    [os.remove(pdbfile) for pdbfile in glob('{0}/*.pdb'.format(ZDOCK_DIR))]


# read ZDOCK_SCORE_FN (if it exists)
if os.path.isfile(ZDOCK_SCORE_FN):
    print("read ZDOCK_SCORE_FN...")
    in_fptr = open(ZDOCK_SCORE_FN)
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
        gene_pdb_dictionary[match.group('gene1')] = string2list(match.group('pdbset1'))
        gene_pdb_dictionary[match.group('gene2')] = string2list(match.group('pdbset2'))
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
        # store rows without Interactome3D results in memory
        fi_no_interactome_set.add((match.group('gene1'),
                                   match.group('swissprot1'),
                                   match.group('gene2'),
                                   match.group('swissprot2')))
# close file
in_fptr.close()

# write ZDOCK_SCORE_FN header: genes, pdbs, top zscore, number of complexes with top score
print("write ZDOCK_SCORE_FN header...")
with open(ZDOCK_SCORE_FN, 'w') as out_fptr:
    out_fptr.write("gene1,pdb_receptor,gene2,pdb_ligand,zdock_score,num_complexes\n")

# for stored FIINT_FN rows
for interaction in fi_no_interactome_set:
    print("for interaction...")
    # for each pair of pdbs (some genes have more than 1) not previously computed
    gene1 = interaction[0]
    swissprot1 = interaction[1]
    gene2 = interaction[2]
    swissprot2 = interaction[3]
    for pdb_tup in itertools.product(gene_pdb_dictionary[gene1], gene_pdb_dictionary[gene2]):
        print("for pdb_tup...")
        # check if .pdb file for each protein exists
        pdb1isfile = os.path.isfile("{0}/{1}.pdb".format(OUT_DIR, pdb_tup[0]))
        pdb2isfile = os.path.isfile("{0}/{1}.pdb".format(OUT_DIR, pdb_tup[1]))

        try:
            if (not pdb1isfile) and (not pdb2isfile):
                # delete any existing .pdb files & download .pdb files for each protein
                [os.remove(pdbfile) for pdbfile in glob('{0}/*.pdb'.format(OUT_DIR))]
                [os.remove(pdbfile) for pdbfile in glob('{0}/*.pdb.gz'.format(OUT_DIR))]
                download_pdb(pdb_tup[0], OUT_DIR)
                download_pdb(pdb_tup[1], OUT_DIR)
            elif pdb1isfile:
                # download .pdb 2
                download_pdb(pdb_tup[1], OUT_DIR)
            elif pdb2isfile:
                # download .pdb 1
                download_pdb(pdb_tup[0], OUT_DIR)
        except ftplib.error_perm as ep:
            print("ERROR RETRIEVING PDB: {0} (SKIPPING)".format(ep))
            break

        # perform zdock execution on .pdb files
        print("perform zdock execution...")
        zdock(pdb_tup[0], pdb_tup[1])

        # record top zdock score and number of rows (complexes) with that score
        print("record top zdock score and number of rows with that score...")
        top_zdock_score = float("-inf")
        num_top_scores = 0
        try:
            in_fptr = open(ZDOCK_OUT_FN)
            while 1:
                line = in_fptr.readline()
                if not line:
                    break
                match = re.match(ZDOCK_SCORE_CAPTURE, line)
                if match and float(match.group('zdock_score')) > float("-inf"):
                    if float(match.group('zdock_score')) >= top_zdock_score:
                        top_zdock_score = float(match.group('zdock_score'))
                        num_top_scores += 1
                    else:
                        break  # zdock scores are listed hi -> lo
            in_fptr.close()
        except FileNotFoundError as fnf:
            print("ERROR CREATING ZDOCK.OUT: {0} (SKIPPING)".format(fnf))
            break

        # remove zdock outfile
        print("remove zdock outfile...")
        os.remove(ZDOCK_OUT_FN)

        # open ZDOCK_SCORE_FN
        print("write to ZDOCK_SCORE_FN...")
        # store genes, pdbs, top zscore, number of complexes with top score
        with open(ZDOCK_SCORE_FN, 'a') as out_fptr:
            out_fptr.write("{0},{1},{2},{3},{4},{5}\n".format(gene1, pdb_tup[0],
                                                              gene2, pdb_tup[1],
                                                              top_zdock_score, num_top_scores))
