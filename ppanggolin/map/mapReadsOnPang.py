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

#local libraries
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.sample import Sample
from ppanggolin.formats import writeFastaGeneFam


def map(pangenome, IDsample, readFile1, readFile2, output, cpu = 1):
    checkPangenomeInfo(pangenome, needFamilies=True, needSamples=True)
    writeFastaGeneFam(pangenome, output, False, "all", False)
    fasta_nuc_pangenome = output+"/all_nucleotide_families.fasta"
    index_files_prefix  = fasta_nuc_pangenome + "_index"
    mapping_sam_file    = output + "/reads_mapping.sam"
    cmd = ["bowtie2-build", "-f", fasta_nuc_pangenome, index_files_prefix]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Index cluster representatives using bowties2...")
    #logging.getLogger().debug(subprocess.run(cmd, capture_output=True))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    
    cmd = ["bowtie2", "--no-unal", "-p", str(cpu),"-k", str(1), "-x", index_files_prefix, "-1", readFile1, "-2", readFile2, "-S", mapping_sam_file]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Mapping reads to cluster representatives using bowties2...")
    #logging.getLogger().debug(subprocess.run(cmd, capture_output=True))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    
    count_map = dict(zip([gf.name for gf in pangenome.geneFamilies], [0] * len(pangenome.geneFamilies)))
    with pysam.AlignmentFile(mapping_sam_file, "r") as aln:
        for read_aln in aln:
            count_map[read_aln.reference_name] += 1
    sample=Sample(IDsample)
    sample.gene_families_map_count = count_map
    pangenome.addSample(sample)

    print("compute : "+str(sample.ID)+"   "+str(sample.ID)+"     "+str(sample.gene_families_map_count["GUT_GENOME218257_CDS_0011"]))

    #change status

def import_count_matrix(pangenome, count_matrix_file, sep = "\t"):
    checkPangenomeInfo(pangenome, needFamilies=True, needSamples=True)
    IDsamples = []
    parse_header= True
    with open(count_matrix_file,"r") as count_matrix:
        for line in count_matrix:
            elements = line.strip("\n").split(sep)
            if parse_header:
                IDsamples = elements[1:]
                initialized_counts = dict(zip([gf.name for gf in pangenome.geneFamilies], [0] * len(pangenome.geneFamilies)))
                for IDsample in IDsamples:
                    sample=Sample(IDsample)
                    sample.gene_families_map_count = initialized_counts.copy()
                    pangenome.addSample(sample)
                    print(IDsample)
                parse_header=False
                continue
            else:
                genefamily =  elements[0]
                genefamily_counts = elements[1:]
                for indexSample, IDsample in enumerate(IDsamples):
                    pangenome.samples[IDsample].gene_families_map_count[genefamily] = float(genefamily_counts[indexSample])

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.reads1 is not None and args.ID is not None:
        map(pangenome = pangenome, IDsample = args.ID, readFile1 = args.reads1, readFile2 = args.reads2, output = args.output, cpu= args.cpu)
    elif args.import_count_matrix is not None:
        import_count_matrix(pangenome = pangenome, count_matrix_file=args.import_count_matrix)
    pangenome.status["samples"]="Computed"
    writePangenome(pangenome, pangenome.file, args.force, show_bar=args.show_prog_bars)

def mapSubparser(subparser):
    parser = subparser.add_parser("map", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Mandatory input/output files", description = "One of the following argument is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")

    mapping = parser.add_argument_group(title = "mapping arguments", description = "")
    mapping.add_argument('--ID', required = False, type = str, help = "sample ID")
    mapping.add_argument('--reads1', required = False, type = str, help = "reads 1 to map against the pangenome gene families")
    mapping.add_argument('--reads2', required = False, type = str, help = "reads 2 to map against the pangenome gene families")

    cout_matrix = parser.add_argument_group(title = "import count matrix", description = "")
    cout_matrix.add_argument('--import_count_matrix', default=None, required = False, type = str, help = "import a already computed count matrix")
    
    #optional = parser.add_argument_group(title = "Optional arguments")
    return parser