#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import time
import argparse
import logging

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_file_name, check_option_workflow
from ppanggolin.annotate.annotate import annotate_pangenome, read_annotations, get_gene_sequences_from_fastas
from ppanggolin.cluster.cluster import clustering, read_clustering
from ppanggolin.graph.makeGraph import compute_neighbors_graph
from ppanggolin.nem.rarefaction import make_rarefaction_curve
from ppanggolin.nem.partition import partition
from ppanggolin.formats.writeBinaries import write_pangenome
from ppanggolin.formats.writeFlat import write_flat_files
from ppanggolin.figures.ucurve import draw_ucurve
from ppanggolin.figures.tile_plot import draw_tile_plot
from ppanggolin.info.info import print_info
from ppanggolin.mod.module import predict_modules


"""a global workflow that does everything in one go."""


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    check_option_workflow(args)
    pangenome = Pangenome()
    filename = mk_file_name(args.basename, args.output, args.force)
    writing_time, anno_time, clust_time, mod_time, desc_time = (None, None, None, None, None)
    if args.anno:  # if the annotations are provided, we read from it
        start_anno = time.time()
        read_annotations(pangenome, args.anno, cpu=args.cpu, disable_bar=args.disable_prog_bar)
        anno_time = time.time() - start_anno
        start_writing = time.time()
        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        writing_time = time.time() - start_writing
        if args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is None:
            raise Exception("The gff/gbff provided did not have any sequence informations, "
                            "you did not provide clusters and you did not provide fasta file. "
                            "Thus, we do not have the information we need to continue the analysis.")

        elif args.clusters is None and pangenome.status["geneSequences"] == "No" and args.fasta is not None:
            get_gene_sequences_from_fastas(pangenome, args.fasta)
        start_clust = time.time()
        if args.clusters is not None:
            read_clustering(pangenome, args.clusters, disable_bar=args.disable_prog_bar)

        elif args.clusters is None:  # we should have the sequences here.
            clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, disable_bar=args.disable_prog_bar)
        clust_time = time.time() - start_clust
    elif args.fasta is not None:
        start_anno = time.time()
        annotate_pangenome(pangenome, args.fasta, args.tmpdir, args.cpu, disable_bar=args.disable_prog_bar)
        anno_time = time.time() - start_anno
        start_writing = time.time()
        write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
        writing_time = time.time() - start_writing
        start_clust = time.time()
        clustering(pangenome, args.tmpdir, args.cpu, defrag=not args.no_defrag, disable_bar=args.disable_prog_bar)
        clust_time = time.time() - start_clust

    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    start_graph = time.time()
    compute_neighbors_graph(pangenome, disable_bar=args.disable_prog_bar)
    graph_time = time.time() - start_graph

    start_part = time.time()
    partition(pangenome, tmpdir=args.tmpdir, kval=args.nb_of_partitions, cpu=args.cpu,
              disable_bar=args.disable_prog_bar)
    part_time = time.time() - start_part

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    start_mods = time.time()
    predict_modules(pangenome=pangenome, tmpdir=args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar)
    mod_time = time.time() - start_mods

    start_writing = time.time()
    write_pangenome(pangenome, filename, args.force, disable_bar=args.disable_prog_bar)
    writing_time = writing_time + time.time() - start_writing

    if args.rarefaction:
        make_rarefaction_curve(pangenome, args.output, args.tmpdir, cpu=args.cpu, disable_bar=args.disable_prog_bar)
    if 1 < len(pangenome.organisms) < 5000:
        draw_tile_plot(pangenome, args.output, nocloud=False if len(pangenome.organisms) < 500 else True)
    draw_ucurve(pangenome, args.output)

    start_desc = time.time()
    write_flat_files(pangenome, args.output, args.cpu, csv=True, gene_pa=True, gexf=True, light_gexf=True,
                     projection=True, stats=True, json=True, partitions=True, modules=True)
    desc_time = time.time() - start_desc

    logging.getLogger().info(f"Annotation took : {round(anno_time, 2)} seconds")
    logging.getLogger().info(f"Clustering took : {round(clust_time, 2)} seconds")
    logging.getLogger().info(f"Building the graph took : {round(graph_time, 2)} seconds")
    logging.getLogger().info(f"Partitioning the pangenome took : {round(part_time, 2)} seconds")
    logging.getLogger().info(f"Predicting modules took : {round(mod_time, 2)} seconds")
    logging.getLogger().info(f"Writing the pangenome data in HDF5 took : {round(writing_time, 2)} seconds")
    logging.getLogger().info(f"Writing descriptive files for the pangenome took : {round(desc_time, 2)} seconds")
    print_info(filename, content=True)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("panmodule", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Input arguments", description="The possible input arguments :")
    required.add_argument('--fasta', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). "
                               "One line per organism. This option can be used alone.")
    required.add_argument('--anno', required=False, type=str,
                          help="A tab-separated file listing the organism names, "
                               "and the gff filepath of its annotations (the gffs can be compressed). "
                               "One line per organism. This option can be used alone "
                               "IF the fasta sequences are in the gff files, otherwise --fasta needs to be used.")
    required.add_argument("--clusters", required=False, type=str,
                          help="a tab-separated file listing the cluster names, the gene IDs, "
                               "and optionally whether they are a fragment or not.")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=str,
                          default="ppanggolin_output" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                      time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory")
    optional.add_argument("--basename", required=False, default="pangenome", help="basename for the output file")
    optional.add_argument("--rarefaction", required=False, action="store_true",
                          help="Use to compute the rarefaction curves (WARNING: can be time consuming)")
    optional.add_argument("-K", "--nb_of_partitions", required=False, default=-1, type=int,
                          help="Number of partitions to use. Must be at least 2. "
                               "If under 2, it will be detected automatically.")
    optional.add_argument("--defrag", required=False, action="store_true", help=argparse.SUPPRESS)
    # This ensures compatibility with workflows built with the old option "defrag" when it was not the default
    optional.add_argument("--no_defrag", required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with their non-fragmented gene family.")

    return parser
