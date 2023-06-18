#!/usr/bin/env python3
# coding:utf-8

# default libraries
import sys
import os
import gzip
import argparse
from io import TextIOWrapper
from pathlib import Path
from typing import TextIO, Union, BinaryIO

import networkx as nx
import pkg_resources
from numpy import repeat
import logging

from scipy.sparse import csc_matrix

from ppanggolin.geneFamily import GeneFamily


def check_log(name: str) -> TextIO:
    """Check if the output log is writable

    :param name: Path to the log output

    :return: output for log
    """
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        try:
            log_file = open(name, "w")
        except IOError:
            raise IOError("The given log file does not appear.")
        except Exception:
            raise Exception("An unexpected error happened with your logfile. Please check if he is accessible."
                            "If everything looks good, please report an issue on our GitHub.")
        else:
            return log_file


def check_tsv_sanity(tsv):
    """ Check if the given tsv is readable for the next PPanGGOLiN step

    :param tsv: Path to the input tsv
    """
    f = open(tsv, "r")
    name_set = set()
    duplicated_names = set()
    non_existing_files = set()
    for line in f:
        elements = [el.strip() for el in line.split("\t")]
        if len(elements) <= 1:
            raise Exception(f"No tabulation separator found in given file: {tsv}")
        if " " in elements[0]:
            raise Exception(f"Your genome names contain spaces (The first encountered genome name that had this string:"
                            f" '{elements[0]}'). To ensure compatibility with all of the dependencies of PPanGGOLiN "
                            f"this is not allowed. Please remove spaces from your genome names.")
        old_len = len(name_set)
        name_set.add(elements[0])
        if len(name_set) == old_len:
            duplicated_names.add(elements[0])
        if not os.path.exists(elements[1]):
            non_existing_files.add(elements[1])
    if len(non_existing_files) != 0:
        raise Exception(f"Some of the given files do not exist. The non-existing files are the following : "
                        f"'{' '.join(non_existing_files)}'")
    if len(duplicated_names) != 0:
        raise Exception(f"Some of your genomes have identical names. The duplicated names are the following : "
                        f"'{' '.join(duplicated_names)}'")


def check_input_files(anno: str = None, pangenome: str = None, fasta: str = None):
    """ Checks if the provided input files exist and are of the proper format

    :param anno: Path to the annotation file
    :param pangenome: Path to the pangenome hdf5 file
    :param fasta: path to the fasta file
    """
    if pangenome is not None and not os.path.exists(pangenome):
        raise FileNotFoundError(f"No such file or directory: '{pangenome}'")

    if anno is not None:
        if not os.path.exists(anno):
            raise FileNotFoundError(f"No such file or directory: '{anno}'")
        check_tsv_sanity(anno)

    if fasta is not None:
        if not os.path.exists(fasta):
            raise FileNotFoundError(f"No such file or directory: '{fasta}'")
        check_tsv_sanity(fasta)


def set_verbosity_level(args):
    """Set the verbosity level

    :param args: argument pass by command line
    """
    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True

        logging.basicConfig(stream=args.log, level=level,
                            format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.getLogger().info("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger().info("PPanGGOLiN version: " + pkg_resources.get_distribution("ppanggolin").version)


def jaccard_similarities(mat: csc_matrix, jaccard_similarity_th) -> csc_matrix:
    """ Compute the jaccard similarities

    :param mat:
    :param jaccard_similarity_th: threshold

    :return:
    """
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat
    # for rows
    aa = repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]
    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)
    similarities.data[similarities.data < jaccard_similarity_th] = 0
    similarities.eliminate_zeros()
    return similarities


def read_compressed_or_not(file_or_file_path: Union[str, BinaryIO, TextIOWrapper, TextIO]) -> Union[TextIOWrapper,
                                                                                                    BinaryIO, TextIO]:
    """
    Reads a file object or file path, uncompresses it, if need be.

    :param file_or_file_path: Path to the input file

    :return: TextIO object in read only
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return file
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        return TextIOWrapper(gzip.open(filename=file, mode="r"))
    else:
        file.close()
        file = open(file.name, "r")
        return file


def write_compressed_or_not(file_path: str, compress: bool = False) -> Union[gzip.GzipFile, TextIO]:
    """
    Create a file-like object, compressed or not.

    :param file_path: Path to the file
    :param compress: Compress the file in .gz

    :return: file-like object, compressed or not
    """
    if compress:
        return gzip.open(file_path + ".gz", mode="wt")
    else:
        return open(file_path, "w")


def is_compressed(file_or_file_path: Union[str, TextIO, gzip.GzipFile]):
    """ Checks is a file, or file path given is compressed or not

    :param file_or_file_path: Input file

    :return: Get if the file is compressed
    """
    file = file_or_file_path
    if isinstance(file, str):
        file = open(file, "rb")
    else:
        try:
            file = open(file.name, "rb")
        except AttributeError:
            return False
    if file.read(2).startswith(b'\x1f\x8b'):
        return True
    file.close()
    return False


def mk_outdir(output, force):
    """ Create a directory at the given output if it doesn't exist already

    :param output: Path where to create directory
    :param force: Force to write in the directory

    :raise FileExistError: The current path already exist and force is false
    """
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


def mk_file_name(basename: str, output: str, force: bool = False) -> Path:
    """Returns a usable filename for a ppanggolin output file, or crashes.

    :param basename: basename for the file
    :param output: Path to save the file
    :param force: Force to write the file

    :return: Path to the file
    """
    filename = Path(output + "/" + basename)
    if filename.suffix != ".h5":
        filename = filename.with_suffix(".h5")

    mk_outdir(output, force)

    if filename.exists() and not force:
        raise FileExistsError(f"{filename.name} already exists. Use -f if you want to overwrite the file")
    return filename


def restricted_float(x) -> float:
    """Decrease the choice possibility of float in argparse

    :param x: given float by user

    :return: given float if it is acceptable

    :raise argparse.ArgumentTypeError: The float is not acceptable
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def min_one(x) -> int:
    """Check if the given int is superior to one

    :param x: given float by user

    :return: given float if it is acceptable

    :raise argparse.ArgumentTypeError: The float is not acceptable
    """
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("%r is inferior to 1" % (x,))
    return x


def connected_components(g: nx.Graph, removed: set, weight: float):
    """
    Yields subgraphs of each connected component you get when filtering edges based on the given weight.

    :param g: Subgraph
    :param removed: removed node
    :param weight: threshold to remove node or not
    """
    for v in g.nodes:
        if v not in removed:
            c = set(_plain_bfs(g, v, removed, weight))
            yield c
            removed.update(c)


def _plain_bfs(g: nx.Graph, source: GeneFamily, removed: set, weight: float):
    """
    A fast BFS node generator, copied from networkx then adapted to the current use case

    :param g: graph with the nodes
    :param source: current node
    :param removed: set of removed nodes
    :param weight:threshold to remove node or not
    """

    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in removed:
                yield v
                removed.add(v)

                for n in g.neighbors(v):
                    if n not in removed:
                        edge_genes_v = g[v][n]["genes"][v]
                        edge_genes_n = g[v][n]["genes"][n]
                        # if the edge is indeed existent for most genes of both families, we use it
                        if len(edge_genes_n) / len(g.nodes[n]["genes"]) >= weight and len(edge_genes_v) / len(
                                g.nodes[v]["genes"]) >= weight:
                            nextlevel.add(n)


def add_gene(obj, gene, fam_split: bool = True):
    """

    :param obj:
    :param gene:
    :param fam_split:
    """
    if fam_split:
        try:
            obj["genes"][gene.family].add(gene)
        except KeyError:
            try:
                obj["genes"][gene.family] = {gene}
            except KeyError:
                obj["genes"] = {gene.family: {gene}}
    else:
        try:
            obj["genes"].add(gene)
        except KeyError:
            obj["genes"] = {gene}


def check_option_workflow(args):
    """
    Check if the given argument to a workflow command is usable

    :param args: list of arguments
    """
    if args.clusters is not None and not any([args.fasta, args.anno]):
        raise Exception("If you give --clusters option, you must give at least --fasta or --anno")

    if not any([args.fasta, args.anno]):
        raise Exception("At least one of --fasta or --anno must be given")

def get_file_size(f):
    old_file_position = f.tell()
    f.seek(0, os.SEEK_END)
    size = f.tell()
    f.seek(old_file_position, os.SEEK_SET)
    return size

reverse_sign = lambda x: '-' if (x=='+') else '+'

def canonical(node_list: str):
    """
    Returns the canonical representation of a node list. 
    We define the canonical as the min between a list eg: 684619+,684620+,684618+ and its reverse eg: 684618-,684620-,684619-
    wrt to the first unsigned value (here 684619 or 684618). In this case we return 684618-,684620-,684619-
    """
    splitted_node_list = node_list.split(',')
    # if only one node, sign is always "+"
    if len(splitted_node_list) == 1:
        return splitted_node_list[0][:-1]+"+"
    else:
        if int(splitted_node_list[0][:-1]) < int(splitted_node_list[-1][:-1]): 
            return node_list
        # else
        rev_list = ','.join([val[:-1]+reverse_sign(val[-1]) for val in reversed(splitted_node_list)])
        return rev_list
        
def canonical_tiny(node_list: list):
    """
    Returns the canonical representation of a node list.
    For smaller lists (reads mapping paths).
    Check the sum of the difference direction of each intervalles between successive nodes
    """
    if len(node_list) == 1 or sum([t - s > 0 for s, t in zip(node_list, node_list[1:])])/(len(node_list)-1) > 0.5:
        return node_list
    return node_list[::-1]

def filter_fastq(input_path, output_path, headers_to_keep):
    with open(output_path, 'w') as output_fastq:
        with read_compressed_or_not(input_path) as input_fastq:
            next_lines_to_print = 0
            for line in input_fastq:
                if next_lines_to_print < 0:
                    next_lines_to_print += 1
                    continue
                elif next_lines_to_print > 0:
                    output_fastq.write(line)
                    next_lines_to_print -= 1
                elif line.startswith("@"):
                    if line.split(" ")[0][1:].strip() in headers_to_keep:
                        output_fastq.write(line)
                        next_lines_to_print = 3
                    else:
                        next_lines_to_print = -3
                else:
                    print("problem " + line)
