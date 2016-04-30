from lxml import html
import requests
import urllib
import os
import codecs
import sys
import re
from graphviz import Graph
from graphviz import Digraph
from collections import defaultdict
from markdown import markdown,markdownFromFile
import pdfkit

PFAM_STRING = 'Pfam'
REPORT_TEMPLATE_DIR = '{0}/templates'.format(os.path.dirname(os.path.realpath(__file__)))
REPORT_TEMPLATE_FN = '{0}/ReportTemplate.html'.format(REPORT_TEMPLATE_DIR)
REPORT_TEMPLATE_VIS_DIR = '{0}/vis'.format(REPORT_TEMPLATE_DIR)
FUSION_INTRCN_RESULTS_FN = '{0}/../data/output/FIInteractFusionEvents.txt'.format(
    os.path.dirname(os.path.realpath(__file__)))
INTERACTOME_DETAIL_CAPTURE = '^F\s(?P<fusion1>[a-zA-Z0-9]+)-(?P<fusion2>[a-z-A-Z0-9]+),I\s(?P<intrcn1>[a-z-A-Z0-9]+)-(?P<intrcn2>[a-z-A-Z0-9]+),(?P<gene1>[a-zA-Z0-9]+),(?P<gene1uni>[a-zA-Z0-9]+),.+,(?P<gene2>[a-zA-Z0-9]+),(?P<gene2uni>[a-zA-Z0-9]+),.+,Interaction.+'
FUSION_INTRCN_ONLY_CAPTURE = '^F\s(?P<fusion>[a-zA-Z0-9]+-[a-z-A-Z0-9]+),I\s(?P<intrcn>[a-z-A-Z0-9]+-[a-z-A-Z0-9]+),.+'
SVG_DIR = '{0}/../data/output/vis/svg'.format(os.path.dirname(os.path.realpath(__file__)))

# Parse Input
f_i_dict = dict()
rep_i_set = set()

if os.path.isfile(FUSION_INTRCN_RESULTS_FN):
    print('reading {0}...'.format(FUSION_INTRCN_RESULTS_FN))
    in_fptr = open(FUSION_INTRCN_RESULTS_FN)
    while 1:
        line = in_fptr.readline()
        if not line:
            break
        match = re.match(FUSION_INTRCN_ONLY_CAPTURE, line)
        if match:
            f = match.group('fusion')
            i = match.group('intrcn')
            print('adding {0} to {1}...'.format(i,f))
            f_i_dict[f] = f_i_dict[f] + [i] if f_i_dict.get(f) is not None else [i]
        match = re.match(INTERACTOME_DETAIL_CAPTURE, line)
        if match:
            print('interaction report generated for {0}'.format(i))
            rep_i_set.add(i)

    # close file
    in_fptr.close()

    # Generate Fusion-Interaction Mapping Image
    dot = Digraph(comment='Fusion-Interaction Mapping',format='svg')

    # Nodes
    for fusion, intrcn_list in f_i_dict.items():
        print('fusion: {0}'.format(fusion))
        print('intrcn_list: {0}'.format(intrcn_list))
        fus = Digraph(name=fusion)
        fus.node(fusion,fusion,style='filled',fillcolor='red')
        for intrcn in intrcn_list:
            print('intrcn: {0}'.format(intrcn))
            if intrcn in rep_i_set:
                fus.node(intrcn,intrcn,shape='box',style='filled',fillcolor='yellow',color='cyan')
            else:
                fus.node(intrcn,intrcn,shape='box',style='filled',fillcolor='gray')
            # Edges
            fus.edge(fusion,intrcn)
        dot.subgraph(fus)

    fusion_intrcn_mapping_svg_path = '{0}/fusion_intrcn_mapping'.format(SVG_DIR)
    dot.render(fusion_intrcn_mapping_svg_path,view=False)
    os.system('ccomps -x {0} | dot | gvpack -array15 | neato -Tsvg -n2 -o {0}_sep.svg'.format(fusion_intrcn_mapping_svg_path))
    fusion_intrcn_mapping_svg_path = '{0}.svg'.format(fusion_intrcn_mapping_svg_path)

    #Note: Use `convert -density 1200 -resize 200x200 source.svg target.png` to redraw as png

    print('fusion-interaction mapping image has been written to {0}.'.format(fusion_intrcn_mapping_svg_path))