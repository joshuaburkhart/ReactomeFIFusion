from lxml import html
import requests
import urllib
import os
import codecs
import re
from graphviz import Digraph
from collections import defaultdict
from markdown import markdown,markdownFromFile
import pdfkit

PFAM_STRING = 'Pfam'
FUSION_INTRCN_RESULTS_FN = '{0}/../data/output/FIInteractFusionEvents.txt'.format(
    os.path.dirname(os.path.realpath(__file__)))
FUSION_INTRCN_CAPTURE = '^F\s(?P<fusion1>[a-zA-Z0-9]+)-(?P<fusion2>[a-z-A-Z0-9]+),I\s(?P<intrcn1>[a-z-A-Z0-9]+)-(?P<intrcn2>[a-z-A-Z0-9]+),(?P<gene1>[a-zA-Z0-9]+),(?P<gene1uni>[a-zA-Z0-9]+),.+,(?P<gene2>[a-zA-Z0-9]+),(?P<gene2uni>[a-zA-Z0-9]+),.+,Interaction.+'

# Parse Input

interactions = list()

if os.path.isfile(FUSION_INTRCN_RESULTS_FN):
    print('reading {0}...'.format(FUSION_INTRCN_RESULTS_FN))
    in_fptr = open(FUSION_INTRCN_RESULTS_FN)
    while 1:
        line = in_fptr.readline()
        if not line:
            break
        match = re.match(FUSION_INTRCN_CAPTURE, line)
        if match:
            intrcn = dict()

            print(line)

            print('fusion1: {0}'.format(match.group('fusion1')))
            print('fusion2: {0}'.format(match.group('fusion2')))
            print('intrcn1: {0}'.format(match.group('intrcn1')))
            print('intrcn2: {0}'.format(match.group('intrcn2')))
            print('gene1: {0}'.format(match.group('gene1')))
            print('gene1uni: {0}'.format(match.group('gene1uni')))
            print('gene2: {0}'.format(match.group('gene2')))
            print('gene2uni: {0}'.format(match.group('gene2uni')))

            intrcn_set = {match.group('intrcn1'), match.group('intrcn2')}
            fusion_set = {match.group('fusion1'), match.group('fusion2')}

            # store map each gene to fusion event
            intrcn['fusion_only_gene'] = fusion_set.difference(intrcn_set).pop() if len(
                fusion_set) > 1 else fusion_set.pop()
            intrcn['intrcn_only_gene'] = intrcn_set.difference(fusion_set).pop() if len(
                intrcn_set) > 1 else intrcn_set.pop()
            intrcn['intrcn_and_fusion_gene'] = {match.group('fusion1'),
                                                match.group('fusion2')} \
                .intersection({match.group('intrcn1'),
                               match.group('intrcn2')}).pop()
            intrcn[match.group('gene1')] = match.group('gene1uni')
            intrcn[match.group('gene2')] = match.group('gene2uni') if match.group('gene1') != match.group(
                'gene2') else [match.group('gene1uni'), match.group('gene2uni')]

            interactions.append(intrcn)
    # close file
    in_fptr.close()

# Loop over Interactions
for intrcn in interactions:
    try:
        fusion_only_gene = intrcn['fusion_only_gene']
        intrcn_only_gene = intrcn['intrcn_only_gene']
        intrcn_and_fusion_gene = intrcn['intrcn_and_fusion_gene']
        intrcn_only_uni = intrcn[intrcn_only_gene] if intrcn_only_gene != intrcn_and_fusion_gene else \
            intrcn[intrcn_only_gene][0]
        intrcn_and_fusion_uni = intrcn[intrcn_and_fusion_gene] if intrcn_only_gene != intrcn_and_fusion_gene else \
            intrcn[intrcn_and_fusion_gene][1]

        uni_gene = {
            intrcn_only_uni: intrcn_only_gene,
            intrcn_and_fusion_uni: intrcn_and_fusion_gene
        }

        # Interactome3D
        interactome3D_url = 'http://interactome3d.irbbarcelona.org/interaction.php?ids={0};{1}&dataset=human&rs=True&connect=1' \
            .format(*(sorted([intrcn_only_uni, intrcn_and_fusion_uni])))

        page = requests.get(interactome3D_url)
        tree = html.fromstring(page.content)

        print(interactome3D_url)

        intrcn_cplx_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[5]/text()'
        cplx1_chain_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[6]/text()'
        cplx1_start_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[1]/div/table/tr/td[1]/text()'
        cplx1_uni_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[1]/div/table/tr/td[2]/text()'
        cplx1_end_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[1]/div/table/tr/td[3]/text()'
        cplx1_from_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[8]/text()'
        cplx1_to_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[9]/text()'
        cplx2_chain_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[12]/text()'
        cplx2_start_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[2]/div/table/tr/td[1]/text()'
        cplx2_uni_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[2]/div/table/tr/td[2]/text()'
        cplx2_end_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[1]/div[2]/div/table/tr/td[3]/text()'
        cplx2_from_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[14]/text()'
        cplx2_to_xpath = '/html/body/div[3]/div[2]/div[2]/div[2]/div[3]/form/table/tr[2]/td[15]/text()'

        intrcn_cplx = tree.xpath(intrcn_cplx_xpath)[0]
        cplx1_chain = tree.xpath(cplx1_chain_xpath)[0]
        cplx1_start = tree.xpath(cplx1_start_xpath)[0]
        cplx1_uni = tree.xpath(cplx1_uni_xpath)[0]
        cplx1_end = tree.xpath(cplx1_end_xpath)[0]
        cplx1_from = tree.xpath(cplx1_from_xpath)[0]
        cplx1_to = tree.xpath(cplx1_to_xpath)[0]
        cplx2_chain = tree.xpath(cplx2_chain_xpath)[0]
        cplx2_start = tree.xpath(cplx2_start_xpath)[0]
        cplx2_uni = tree.xpath(cplx2_uni_xpath)[0]
        cplx2_end = tree.xpath(cplx2_end_xpath)[0]
        cplx2_from = tree.xpath(cplx2_from_xpath)[0]
        cplx2_to = tree.xpath(cplx2_to_xpath)[0]

        # Create Directories
        PDB_DIR = '{0}/../data/output/reports/{1}/src/pdb'.format(os.path.dirname(os.path.realpath(__file__)),
                                                                  intrcn_cplx)
        PNG_DIR = '{0}/../data/output/reports/{1}/src/png'.format(os.path.dirname(os.path.realpath(__file__)),
                                                                  intrcn_cplx)
        MD_DIR = '{0}/../data/output/reports/{1}/md'.format(os.path.dirname(os.path.realpath(__file__)), intrcn_cplx)
        HTML_DIR = '{0}/../data/output/reports/{1}/html'.format(os.path.dirname(os.path.realpath(__file__)),
                                                                intrcn_cplx)
        PDF_DIR = '{0}/../data/output/reports/{1}/pdf'.format(os.path.dirname(os.path.realpath(__file__)), intrcn_cplx)

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
        pdb_gz_path = '{0}/{1}.pdb.gz'.format(PDB_DIR, intrcn_cplx)
        if os.path.exists(pdb_gz_path):
            os.remove(pdb_gz_path)
        urllib.request.urlretrieve('http://files.rcsb.org/download/{0}.pdb.gz'.format(intrcn_cplx.upper()), pdb_gz_path)
        pdb_path = '{0}/{1}.pdb'.format(PDB_DIR, intrcn_cplx)
        if os.path.exists(pdb_path):
            os.remove(pdb_path)
        os.system('gunzip {0}'.format(pdb_gz_path))

        # Generate Jmol Image
        jmol_png_path = '{0}/{1}.png'.format(PNG_DIR, intrcn_cplx)
        if os.path.exists(jmol_png_path):
            os.remove(jmol_png_path)
        os.system('java -Djava.awt.headless=false -Xmx512m -jar "/Applications/Jmol.jar" -onj \
    "load file {0}; \
    select all; \
    color grey; \
    select:{1}; \
    color yellow; \
    select:{2}; \
    color aqua; \
    show BEST ROTATION; \
    write {3}"'.format(pdb_path, cplx1_chain, cplx2_chain, jmol_png_path))

        # PDBe
        pdbe_root_url = 'http://www.ebi.ac.uk/pdbe/entry/pdb/{0}'.format(intrcn_cplx)

        # Biology
        page = requests.get('{0}/{1}'.format(pdbe_root_url, 'biology'))
        tree = html.fromstring(page.content)
        rctn_xpath = '/html/body/div[2]/div/section/div[3]/div/section[1]/div/div/section[1]/div/text()'
        rctn = tree.xpath(rctn_xpath)[0].strip() if len(
            tree.xpath(rctn_xpath)) > 0 else 'Biological process not assigned.'

        # Experiment
        page = requests.get('{0}/{1}'.format(pdbe_root_url, 'experiment'))
        tree = html.fromstring(page.content)
        dscrp_xpath = '/html/body/div[2]/div/section/div[1]/div[3]/div/section[1]/div/div[1]/div/text()'
        dscrp = tree.xpath(dscrp_xpath)[0].strip() if len(tree.xpath(dscrp_xpath)) > 0 else 'Description not available.'

        # Citation
        page = requests.get('{0}/{1}'.format(pdbe_root_url, 'citations'))
        tree = html.fromstring(page.content)
        cit_xpath = '/html/body/div[2]/div/section/div[1]/div[2]/div[2]/div/div[3]/div[2]/a/@href'
        cit = tree.xpath(cit_xpath)[0].strip() if len(tree.xpath(cit_xpath)) > 0 else 'Citation not listed.'

        # Pfam
        pfam_url = 'http://pfam.xfam.org/protein/{0}'.format(intrcn_and_fusion_uni)
        page = requests.get(pfam_url)
        tree = html.fromstring(page.content)

        # Domain Coordinates
        dmn_table_rows_xpath = '/html/body/div[5]/div[5]/div[1]/div[2]/div[1]/div[2]/table[2]//tr[@class="odd"]'
        dmn_table_rows = tree.xpath(dmn_table_rows_xpath)

        dmns = list()

        for tr in dmn_table_rows:
            if tr[0].text == PFAM_STRING:
                dmns.append({
                    'domain': tr[1][0].text,
                    'start': tr[2].text,
                    'end': tr[3].text
                })

        # COSMIC

        cosmic_fusion_url = 'http://cancer.sanger.ac.uk/cosmic/stats/fusion'
        page = requests.get(cosmic_fusion_url)
        tree = html.fromstring(page.content)

        fusion_link_xpath = './/a[contains(text(),"{0}")][contains(text(),"{1}")]/@href'.format(fusion_only_gene,
                                                                                                intrcn_and_fusion_gene)

        print(fusion_only_gene)
        print(intrcn_and_fusion_gene)

        fusion_link = tree.xpath(fusion_link_xpath)[0]

        page = requests.get(fusion_link)
        tree = html.fromstring(page.content)

        mutation_id_xpath = '/html/body/div/div[3]/div/div[1]/div[2]/div[3]/table/tr/td[1]/a/text()'
        gene1_xpath = '/html/body/div/div[3]/div/div[1]/div[2]/div[3]/table/tr/td[2]/a/text()'
        gene1brk_xpath = '/html/body/div/div[3]/div/div[1]/div[2]/div[3]/table/tr/td[4]/text()'
        gene2_xpath = '/html/body/div/div[3]/div/div[1]/div[2]/div[3]/table/tr/td[6]/a/text()'
        gene2brk_xpath = '/html/body/div/div[3]/div/div[1]/div[2]/div[3]/table/tr/td[8]/text()'

        mutation_id = tree.xpath(mutation_id_xpath)[0]
        cosmic_gene1 = tree.xpath(gene1_xpath)[0]
        gene1brk = tree.xpath(gene1brk_xpath)[0]
        cosmic_gene2 = tree.xpath(gene2_xpath)[0]
        gene2brk = tree.xpath(gene2brk_xpath)[0]

        gene_brk = {
            cosmic_gene1: gene1brk,
            cosmic_gene2: gene2brk
        }

        aa_brk = {
            cosmic_gene1: int(int(gene_brk[cosmic_gene1].split('+')[-1:][0].split('-')[-1:][0]) / 3),
            cosmic_gene2: int(int(gene_brk[cosmic_gene2].split('+')[-1:][0].split('-')[-1:][0]) / 3)
        }

        # Generate Interaction Schematic Image
        dot = Digraph(comment='Interaction Schematic')

        # Nodes
        dot.node('A', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], cplx1_start, int(cplx1_from) - 1), shape='box')
        dot.node('B', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], cplx1_from, cplx1_to), color='yellow', shape='box')
        dot.node('C', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], int(cplx1_to) + 1, cplx1_end), shape='box')
        dot.node('E', '{0} Chain {1}'.format(intrcn_cplx, cplx1_chain), color='yellow', shape='box')
        dot.node('F', '{0} Chain {1}'.format(intrcn_cplx, cplx2_chain), color='cyan', shape='box')
        dot.node('G', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], cplx2_start, int(cplx2_from) - 1), shape='box')
        dot.node('H', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], cplx2_from, cplx2_to), color='cyan', shape='box')
        dot.node('I', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], int(cplx2_to) + 1, cplx2_end), shape='box')

        # Edges
        dot.edges(['AB', 'BC', 'GH', 'HI', 'EF'])
        dot.edge('B', 'E', constraint='false', color='gray')
        dot.edge('H', 'F', constraint='false', color='gray')

        dot.format = 'png'
        intrcn_schematic_png_path = '{0}/{1}_scheme'.format(PNG_DIR, intrcn_cplx)
        dot.render(intrcn_schematic_png_path, view=False)
        intrcn_schematic_png_path = '{0}.png'.format(intrcn_schematic_png_path)

        # Generate Fusion Effect Schematic Image
        style = {}
        style = defaultdict(lambda: '',style)

        fillcolor = {}
        fillcolor = defaultdict(lambda: 'lightgrey',fillcolor)
        # Interaction Complex 1
        if (uni_gene[cplx1_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx1_from)-1) or \
                (uni_gene[cplx1_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx1_from)-1):
            # Node A
            style['A'] = 'filled'
            fillcolor['A'] = 'red'
            # Node B
            style['B'] = 'filled'
            fillcolor['B'] = 'pink'
            # Node C
            style['C'] = 'filled'
            fillcolor['C'] = 'pink'
            # Node E (interaction chain)
            style['E'] = 'filled'
            fillcolor['E'] = 'pink'
        elif (uni_gene[cplx1_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx1_to)) or \
                (uni_gene[cplx1_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx1_to)):
            # Node B
            style['B'] = 'filled'
            fillcolor['B'] = 'red'
            # Node C
            style['C'] = 'filled'
            fillcolor['C'] = 'pink'
            # Node E (interaction chain)
            style['E'] = 'filled'
            fillcolor['E'] = 'pink'
        elif (uni_gene[cplx1_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx1_end)) or \
                (uni_gene[cplx1_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx1_end)):
            # Node C
            style['C'] = 'filled'
            fillcolor['C'] = 'red'
        # Interaction Complex 2
        elif (uni_gene[cplx2_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx2_from)-1) or \
                (uni_gene[cplx2_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx2_from)-1):
            # Node G
            style['G'] = 'filled'
            fillcolor['G'] = 'red'
            # Node H
            style['H'] = 'filled'
            fillcolor['H'] = 'pink'
            # Node I
            style['I'] = 'filled'
            fillcolor['I'] = 'pink'
            # Node F (interaction chain)
            style['F'] = 'filled'
            fillcolor['F'] = 'pink'
        elif (uni_gene[cplx2_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx2_to)) or \
                (uni_gene[cplx2_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx2_to)):
            # Node H
            style['H'] = 'filled'
            fillcolor['H'] = 'red'
            # Node I
            style['I'] = 'filled'
            fillcolor['I'] = 'pink'
            # Node F (interaction chain)
            style['F'] = 'filled'
            fillcolor['F'] = 'pink'
        elif (uni_gene[cplx2_uni] == cosmic_gene1 and aa_brk[cosmic_gene1] < int(cplx2_end)) or \
                (uni_gene[cplx2_uni] == cosmic_gene2 and aa_brk[cosmic_gene2] < int(cplx2_end)):
            # Node I
            style['I'] = 'filled'
            fillcolor['I'] = 'red'

        dot = Digraph(comment='Fusion Effect Schematic')

        # Nodes
        dot.node('A', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], cplx1_start, int(cplx1_from) - 1), shape='box',style=style['A'],fillcolor=fillcolor['A'])
        dot.node('B', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], cplx1_from, cplx1_to), color='yellow', shape='box',style=style['B'],fillcolor=fillcolor['B'])
        dot.node('C', '{0}: {1}-{2}'.format(uni_gene[cplx1_uni], int(cplx1_to) + 1, cplx1_end), shape='box',style=style['C'],fillcolor=fillcolor['C'])
        dot.node('E', '{0} Chain {1}'.format(intrcn_cplx, cplx1_chain), color='yellow', shape='box',style=style['E'],fillcolor=fillcolor['E'])
        dot.node('F', '{0} Chain {1}'.format(intrcn_cplx, cplx2_chain), color='cyan', shape='box',style=style['F'],fillcolor=fillcolor['F'])
        dot.node('G', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], cplx2_start, int(cplx2_from) - 1), shape='box',style=style['G'],fillcolor=fillcolor['G'])
        dot.node('H', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], cplx2_from, cplx2_to), color='cyan', shape='box',style=style['H'],fillcolor=fillcolor['H'])
        dot.node('I', '{0}: {1}-{2}'.format(uni_gene[cplx2_uni], int(cplx2_to) + 1, cplx2_end), shape='box',style=style['I'],fillcolor=fillcolor['I'])

        # Edges
        dot.edges(['AB', 'BC', 'GH', 'HI', 'EF'])
        dot.edge('B', 'E', constraint='false', color='gray')
        dot.edge('H', 'F', constraint='false', color='gray')

        dot.format = 'png'
        fusion_schematic_png_path = '{0}/{1}_fusion_scheme'.format(PNG_DIR, intrcn_cplx)
        dot.render(fusion_schematic_png_path, view=False)
        fusion_schematic_png_path = '{0}.png'.format(fusion_schematic_png_path)

        # Report Content Logic
        fusion_effect = '(No direct effect on interaction {0})'.format(intrcn_cplx)
        if intrcn_and_fusion_gene == cosmic_gene1:
            if intrcn_and_fusion_gene == uni_gene[cplx1_uni]:
                if cplx1_start < gene_brk[cosmic_gene1]:
                    fusion_effect = 'Chain {0} in leading part of fusion product'.format(cplx1_chain)
            elif intrcn_and_fusion_gene == uni_gene[cplx2_uni]:
                if cplx2_start < gene_brk[cosmic_gene1]:
                    fusion_effect = 'Chain {0} in leading part of fusion product'.format(cplx2_chain)
        elif intrcn_and_fusion_gene == cosmic_gene2:
            if intrcn_and_fusion_gene == uni_gene[cplx1_uni]:
                if cplx1_start > gene_brk[cosmic_gene2]:
                    fusion_effect = 'Chain {0} in trailing part of fusion product'.format(cplx1_chain)
            if intrcn_and_fusion_gene == uni_gene[cplx2_uni]:
                if cplx2_start > gene_brk[cosmic_gene2]:
                    fusion_effect = 'Chain {0} in trailing part of fusion product'.format(cplx2_chain)

        domain_string = ""
        for dmn in dmns:
            domain_string += '<li>{0}: {1}-{2}</li>\n'.format(dmn['domain'], dmn['start'], dmn['end'])

        cplx1_gene = uni_gene[cplx1_uni]
        cplx2_gene = uni_gene[cplx2_uni]
        gene1_brk = gene_brk[cosmic_gene1]
        gene2_brk = gene_brk[cosmic_gene2]
        aa1_brk = aa_brk[cosmic_gene1]
        aa2_brk = aa_brk[cosmic_gene2]

        # Write Markdown File
        md_path = '{0}/{1}.md'.format(MD_DIR, intrcn_cplx)
        if os.path.exists(md_path):
            os.remove(md_path)
        with open(md_path, 'a') as out_fptr:
            out_fptr.write(
'''\
<h1 style="text-align:center;">Fusion {fusion_only_gene!s}-{intrcn_and_fusion_gene!s} and Interaction {intrcn_cplx!s}: {fusion_effect!s}</h1>
<table style="border: 0px;">
    <tr style="border: 0px;">
        <td style="width:50%;border: 0px;">
            <h2 style="text-align:center;">Interaction {intrcn_cplx!s}</h2>
        </th>
        <td style="width:50%;border:0px;">
            <h2 style="text-align:center;">Fusion Effect</h2>
        </th>
    </tr>
    <tr style="border: 0px;width:50%;">
        <td style="border: 0px;background: white;">
            <img style="vertical-align:bottom;" src="{intrcn_schematic_png_path!s}"/>
        </td>
        <td style="border: 0px;background: white;">
            <img src="{fusion_schematic_png_path!s}"/>
        </td>
    </tr>
    <tr style="border: 0px;">
        <td style="background:white;border: 0px;vertical-align:top;width:50%">
            <a href="{interactome3D_url!s}">
                <h2 style="text-align:center;">Interactome3D</h2>
            </a>
            Interaction Complex: {intrcn_cplx!s}
            <ul>
                <li>{cplx1_gene!s}
                    <ul>
                        <li>Chain: {cplx1_chain!s}</li>
                        <li>Length: {cplx1_start!s} - {cplx1_end!s}</li>
                        <li>Included in Interaction: {cplx1_from!s} - {cplx1_to!s}</li>
                    </ul>
                </li>
                <li>{cplx2_gene!s}
                    <ul>
                        <li>Chain: {cplx2_chain!s}</li>
                        <li>Length: {cplx2_start!s} - {cplx2_end!s}</li>
                        <li>Included in Interaction: {cplx2_from!s} - {cplx2_to!s}</li>
                    </ul>
                </li>
            </ul>
        </td>
        <td style="background:white;border: 0px;">
            <img src="{jmol_png_path!s}"/>
        </td>
    </tr>
    <tr style="border:0px;">
        <td style="background:white;border:0px;vertical-align:top;">
            <a href="{pfam_url!s}">
                <h2 style="text-align:center;">Pfam</h2>
            </a>
            <ul>
                {domain_string!s}
            </ul>
        </td>
        <td style="background:white;border:0px;vertical-align:top;">
            <a href="{fusion_link!s}">
                <h2 style="text-align:center;">COSMIC</h2>
            </a>
            Mutation: {mutation_id!s}
            <ul>
                <li>First Gene in Fusion: {cosmic_gene1!s}
                    <ul>
                        <li>genomic breakpoint: {gene1_brk!s}</li>
                        <li>AA breakpoint: {aa1_brk!s}</li>
                    </ul>
                </li>
                <li>Second Gene in Fusion: {cosmic_gene2!s}
                    <ul>
                        <li>genomic breakpoint: {gene2_brk!s}</li>
                        <li>AA breakpoint: {aa2_brk!s}</li>
                    </ul>
                </li>
            </ul>
        </td>
    </tr>
    <tr style="border:0px;">
        <td style="vertical-align:top;background:white;border:0px;">
            <a href="{pdbe_root_url!s}">
                <h2 style="text-align:center;">PDBe</h2>
            </a>
            <li>Reaction: {rctn!s}</li>
            <li>Author Description: {dscrp!s}</li>
            <li>Citation: <a href="{cit!s}">{cit!s}</a>
            </li>
        </td>
        <td style="text-align:right;vertical-align:bottom;border:0px;background:white;">
            ReactomeFIFusion Generated Report
        </td>
    </tr>
</table>\
'''.format(**locals()))

        # Generate PDF
        with open(md_path, 'r') as in_fptr:
            html_text = markdown(in_fptr.read(), output_format='html4')

        pdfkit.from_string(html_text, '{0}/{1}.pdf'.format(PDF_DIR,intrcn_cplx))

        # Generate HTML
        markdownFromFile(md_path,'{0}/{1}.html'.format(HTML_DIR,intrcn_cplx))

    except IndexError as ie:
        err_msg = 'ERROR: "IndexError" exception when attempting to generate a report for {0}: {1}'.format(intrcn, ie)
        print(err_msg)

        # Attempt to write Markdown File
        md_path = '{0}/{1}.md'.format(MD_DIR, intrcn_cplx)
        if os.path.exists(md_path):
            os.remove(md_path)
        with open(md_path, 'a') as out_fptr:
            out_fptr.write(err_msg)
