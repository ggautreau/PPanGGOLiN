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

import scipy
import statsmodels

#local libraries
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.sample import Sample
from ppanggolin.compare import Dataset, Comparison

def compare(pangenome, dataset1, dataset2, samples_dataset1, samples_dataset2, paired = False, correct_p_values = True):
    if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    checkPangenomeInfo(pangenome, needFamilies=True, needSamples=True, needComparisons=True)
    print(pangenome.samples)
    print(samples_dataset1)
    print(dataset1)
    print([pangenome.samples[s] for s in samples_dataset1])
    print(samples_dataset2)
    print(dataset2)
    print([pangenome.samples[s] for s in samples_dataset2])
    ds1 = Dataset(dataset1, [pangenome.samples[s] for s in samples_dataset1])
    ds2 = Dataset(dataset2, [pangenome.samples[s] for s in samples_dataset2])
    pangenome.addDataset(ds1)
    pangenome.addDataset(ds2)
    pangenome.addComparison(Comparison(dataset1 + "_versus_" + dataset2, ds1, ds2, paired, correct_p_values))

    pangenome.status["comparisons"] = "Computed"

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.dataset1 is not None and args.dataset2 is not None and args.samples_dataset1 is not None and args.samples_dataset1 is not None:
        compare(pangenome = pangenome, dataset1 = args.dataset1, dataset2 = args.dataset2, samples_dataset1 = args.samples_dataset1, samples_dataset2 = args.samples_dataset2, paired = args.paired, correct_p_values = args.correct_p_values)
        writePangenome(pangenome, pangenome.file, args.force, show_bar = args.show_prog_bars)

def compareSubparser(subparser):
    parser = subparser.add_parser("compare", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Required arguments", description = "All of the following arguments are required :")
    onereq = parser.add_argument_group(title = "Input file", description = "One of the following argument is required :")
    onereq.add_argument('--dataset1', required = True, type = str, help = "name of the first dataset")
    onereq.add_argument('--dataset2', required = True, type = str, help = "name of the second dataset")
    onereq.add_argument('--samples_dataset1', required = True, nargs = "+", type = str, help = "list samples of the first dataset")
    onereq.add_argument('--samples_dataset2', required = True, nargs = "+", type = str, help = "list samples of the second dataset")
    onereq.add_argument('--paired', required = False, type = str, help = "paired datasets")
    onereq.add_argument('--correct_p_values', required = False, type = str, help = "p_values correction")

    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")
    return parser