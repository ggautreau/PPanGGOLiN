#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict
import pysam
import sys
import pdb

#local libraries
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import writeFastaGeneFam

def import_import_eggnog(pangenome, eggnog_file_path):
    checkPangenomeInfo(pangenome, needFamilies=True)
    with open(eggnog_file_path, "r") as eggnog_file:
        for line in eggnog_file:
            if not line.startswith("#"):
                elements = line.strip("\n").split("\t")
                print(elements)
                pangenome.getGeneFamily(elements[0]).short_name      = elements[4]
                pangenome.getGeneFamily(elements[0]).product_name    = elements[12]
                pangenome.getGeneFamily(elements[0]).eggNOG_ortholog = elements[9]
                pangenome.getGeneFamily(elements[0]).COG             = elements[11]
                pangenome.getGeneFamily(elements[0]).GO_terms        = elements[5]
                pangenome.getGeneFamily(elements[0]).KEGG_KOs        = elements[6]
                pangenome.getGeneFamily(elements[0]).BiGG_reactions  = elements[7]
                
                print(pangenome.getGeneFamily(elements[0]).short_name)
                print(pangenome.getGeneFamily(elements[0]).product_name)
                print(pangenome.getGeneFamily(elements[0]).eggNOG_ortholog)
                print(pangenome.getGeneFamily(elements[0]).COG)
                print(pangenome.getGeneFamily(elements[0]).GO_terms)
                print(pangenome.getGeneFamily(elements[0]).KEGG_KOs)
                print(pangenome.getGeneFamily(elements[0]).BiGG_reactions)

def import_kofam(pangenome, kofam_file_path):
    checkPangenomeInfo(pangenome, needFamilies=True)
    with open(kofam_file_path, "r") as kofam_file:
        for line in kofam_file:
            if line.startswith("*"):
                elements = line[2:].strip().split()
                pangenome.getGeneFamily(elements[0]).KEGG_KOs     = elements[1]
                pangenome.getGeneFamily(elements[0]).product_name = " ".join(elements[5:])
                print(pangenome.getGeneFamily(elements[0]).KEGG_KOs)
                print(pangenome.getGeneFamily(elements[0]).product_name)

def launch(args):
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.import_eggnog is not None:
        import_import_eggnog(pangenome = pangenome, eggnog_file_path = args.import_eggnog)
    elif args.import_kofamkoala is not None:
        import_kofam(pangenome = pangenome, kofam_file_path = args.import_kofamkoala)
    pangenome.status["fonctionImported"] = "Computed"
    writePangenome(pangenome, pangenome.file, args.force, show_bar=args.show_prog_bars)

def functionSubparser(subparser):
    parser = subparser.add_parser("function", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Mandatory input/output files", description = "One of the following argument is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")

    cout_matrix = parser.add_argument_group(title = "import functional annotation eggnog", description = "")
    cout_matrix.add_argument('--import_eggnog', default=None, required = False, type = str, help = "import a already computed eggnog output")
    
    cout_matrix = parser.add_argument_group(title = "import function annotation kofam", description = "")
    cout_matrix.add_argument('--import_kofamkoala', default=None, required = False, type = str, help = "import a already computed kofam output")
    
    return parser