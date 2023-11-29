#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import tempfile
import time
import logging
import os
from typing import List, Dict, Tuple, Iterable, Hashable, Iterator, Set
from itertools import chain
from collections import defaultdict
from pathlib import Path

# installed libraries
from tqdm import tqdm
import networkx as nx
import pandas as pd

# local libraries
from ppanggolin.formats import check_pangenome_info
from ppanggolin.genome import Gene, Contig
from ppanggolin.utils import mk_outdir, restricted_float, create_tmpdir, read_compressed_or_not, extract_contig_window
from ppanggolin.pangenome import Pangenome
from ppanggolin.align.alignOnPang import project_and_write_partition, get_input_seq_to_family_with_rep, \
    get_input_seq_to_family_with_all, get_seq_ids
from ppanggolin.region import GeneContext
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.projection.projection import write_gene_to_gene_family


def search_gene_context_in_pangenome(pangenome: Pangenome, output: Path, tmpdir: Path, sequence_file: Path = None,
                                     families: Path = None, transitive: int = 4, identity: float = 0.5,
                                     coverage: float = 0.8, use_representatives: bool = False,
                                     jaccard_threshold: float = 0.85,
                                     window_size: int = 1, no_defrag: bool = False,
                                     cpu: int = 1, graph_format: str = "graphml", disable_bar=True,
                                     translation_table: int = 11, keep_tmp: bool = False):
    """
    Main function to search common gene contexts between sequence set and pangenome families

    :param pangenome: Pangenome containing GeneFamilies to align with sequence set
    :param sequence_file: Path to file containing the sequences
    :param families: Path to file containing families name
    :param output: Path to output directory
    :param tmpdir: Path to temporary directory
    :param transitive: number of genes to check on both sides of a family aligned with an input sequence
    :param identity: minimum identity threshold between sequences and gene families for the alignment
    :param coverage: minimum coverage threshold between sequences and gene families for the alignment
    :param use_representatives: Use representative sequences of families rather than all sequences to align input genes
    :param jaccard_threshold: Jaccard index threshold to filter edges in graph
    :param window_size: Number of genes to consider in the gene context.
    :param no_defrag: do not use the defrag workflow if true
    :param cpu: Number of core used to process
    :param graph_format: Write format of the context graph. Can be graphml or gexf
    :param disable_bar: Allow preventing bar progress print
    :param translation_table: The translation table to use when the input sequences are nucleotide sequences.
    :param keep_tmp: If True, keep temporary files.
    """
    # check statuses and load info
    if sequence_file is not None and pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)

    families_of_interest = set()
    family_2_input_seqid = {}
    if sequence_file is not None:
        # Alignment of sequences on pangenome families
        with read_compressed_or_not(sequence_file) as seqFileObj:
            seq_set, is_nucleotide = get_seq_ids(seqFileObj)

        logging.debug(f"Input sequences are {'nucleotide' if is_nucleotide else 'protein'} sequences")

        with create_tmpdir(main_dir=tmpdir, basename="align_input_seq_tmpdir", keep_tmp=keep_tmp) as new_tmpdir:

            if use_representatives:
                _, seqid2fam = get_input_seq_to_family_with_rep(pangenome, [sequence_file], output,
                                                                new_tmpdir, is_input_seq_nt=is_nucleotide,
                                                                cpu=cpu, no_defrag=no_defrag,
                                                                identity=identity, coverage=coverage,
                                                                translation_table=translation_table,
                                                                disable_bar=disable_bar)
            else:
                _, seqid2fam = get_input_seq_to_family_with_all(pangenome=pangenome,
                                                                sequence_files=[sequence_file],
                                                                output=output, tmpdir=new_tmpdir,
                                                                is_input_seq_nt=is_nucleotide,
                                                                cpu=cpu, no_defrag=no_defrag,
                                                                identity=identity, coverage=coverage,
                                                                translation_table=translation_table,
                                                                disable_bar=disable_bar)

        project_and_write_partition(seqid2fam, seq_set, output)
        write_gene_to_gene_family(seqid2fam, seq_set, output)

        family_2_input_seqid = defaultdict(set)
        for seqid, gf in seqid2fam.items():
            family_2_input_seqid[gf].add(seqid)

        for pan_family in seqid2fam.values():
            families_of_interest.add(pan_family)

    if families is not None:
        with read_compressed_or_not(families) as f:
            for fam_name in f.read().splitlines():
                families_of_interest.add(pangenome.get_gene_family(fam_name))

    # Compute the graph with transitive closure size provided as parameter
    start_time = time.time()

    logging.getLogger().info("Building the graph...")

    gene_context_graph = compute_gene_context_graph(families=families_of_interest, transitive=transitive,
                                                    window_size=window_size, disable_bar=disable_bar)

    logging.getLogger().info(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")

    logging.getLogger().debug(
        f"Context graph made of {nx.number_of_nodes(gene_context_graph)} nodes and "
        f"{nx.number_of_edges(gene_context_graph)} edges")

    compute_edge_metrics(gene_context_graph, jaccard_threshold)

    # Filter graph 
    filter_flag = f'is_jaccard_gene_>_{jaccard_threshold}'

    edges_to_remove = [(n, v) for n, v, d in gene_context_graph.edges(data=True) if not d[filter_flag]]
    gene_context_graph.remove_edges_from(edges_to_remove)

    logging.getLogger().debug(f"Filtering context graph on {filter_flag}")
    logging.getLogger().debug(
        f"Context graph made of {nx.number_of_nodes(gene_context_graph)} nodes and "
        f"{nx.number_of_edges(gene_context_graph)} edges")

    gene_contexts = get_gene_contexts(gene_context_graph, families_of_interest)

    gene_context_graph = make_graph_writable(gene_context_graph)
    out_graph_file = write_graph(gene_context_graph, output, graph_format)

    if len(gene_contexts) != 0:
        logging.getLogger().info(
            f"There are {sum((len(gc) for gc in gene_contexts))} families among {len(gene_contexts)} gene contexts")

        output_file = output / "gene_contexts.tsv"
        export_context_to_dataframe(gene_contexts, family_2_input_seqid, families_of_interest, output_file)

    else:
        logging.getLogger("PPanGGOLiN").info("No gene contexts were found")

    logging.getLogger("PPanGGOLiN").info(f"Computing gene contexts took {round(time.time() - start_time, 2)} seconds")

    return gene_context_graph, out_graph_file


def get_gene_contexts(context_graph: nx.Graph, families_of_interest: Set[GeneFamily]) -> Set[GeneContext]:
    """
    Extract gene contexts from a context graph based on the provided set of gene families of interest.
    
    Gene contexts are extracted from a context graph by identifying connected components. 
    The function filters the connected components based on the following criteria:
    - Remove singleton families (components with only one gene family).
    - Remove components that do not contain any gene families of interest.

    For each remaining connected component, a GeneContext object is created.

    :param context_graph: The context graph from which to extract gene contexts.
    :param families_of_interest: Set of gene families of interest.
    :return: Set of GeneContext objects representing the extracted gene contexts.
    """

    connected_components = nx.connected_components(context_graph)

    # Connected component graph Filtering

    # remove singleton famillies
    connected_components = (component for component in connected_components if len(component) > 1)

    # remove component made only of famillies not initially requested
    connected_components = (component for component in connected_components if component & families_of_interest)

    gene_contexts = set()
    families_in_context = set()

    for i, component in enumerate(connected_components):
        families_in_context |= component
        family_of_interest_of_gc = component & families_of_interest
        gene_context = GeneContext(gc_id=i, families=component, families_of_interest=family_of_interest_of_gc)

        # add gc id to node attribute
        node_attributes = {n: {"gene_context_id": i, "families_of_interest": n in families_of_interest} for n in
                           component}
        nx.set_node_attributes(context_graph, node_attributes)

        gene_contexts.add(gene_context)

    node_not_in_context = set(context_graph.nodes()) - families_in_context
    context_graph.remove_nodes_from(node_not_in_context)

    return gene_contexts


def make_graph_writable(context_graph):
    """
    The original context graph contains ppanggolin objects as nodes and lists and dictionaries in edge attributes.
    Since these objects cannot be written to the output graph,
    this function creates a new graph that contains only writable objects.

    :param context_graph: List of gene context. it includes graph of the context
    """
    def filter_attribute(data: dict):
        """
        Helper function to filter the edge attributes.

        :param data: The edge attribute data.
        :return: A filtered dictionary containing only non-collection attributes.
        """
        return {k: v for k, v in data.items() if type(v) not in [set, dict, list]}

    G = nx.Graph()

    G.add_edges_from((f1.name, f2.name, filter_attribute(d)) for f1, f2, d in context_graph.edges(data=True))

    # convert transitivity dict to str
    edges_with_transitivity_str = {(f1.name, f2.name): str(d['transitivity']) for f1, f2, d in
                                   context_graph.edges(data=True)}

    nx.set_edge_attributes(G, edges_with_transitivity_str, name="transitivity")

    nodes_attributes_filtered = {f.name: filter_attribute(d) for f, d in context_graph.nodes(data=True)}

    # on top of attributes already contained in node of context graph
    # add organisms and genes count that have the family, the partition and if the family was in initially requested 
    nodes_family_data = {f.name: {"organisms": f.number_of_organisms,
                                  "partition": f.named_partition,
                                  "genes": f.number_of_genes} for f in context_graph.nodes()}

    for f, d in G.nodes(data=True):
        d.update(nodes_family_data[f])
        d.update(nodes_attributes_filtered[f])

    return G


def write_graph(G: nx.Graph, output_dir: Path, graph_format: str):
    """
    Write a graph to file in the GraphML format or/and in GEXF format. 

    :param G: Graph to write
    :param output_dir: The output directory where the graph file will be written.
    :param graph_format: Formats of the output graph. Can be graphml or gexf 

    """

    if "graphml" == graph_format:
        out_file = output_dir / "graph_context.graphml"
        logging.info(f'Writting context graph in {out_file}')
        nx.write_graphml_lxml(G, out_file)

    elif "gexf" == graph_format:
        out_file = output_dir / "graph_context.gexf"
        logging.info(f'Writting context graph in {out_file}')
        nx.readwrite.gexf.write_gexf(G, out_file)
    else:
        raise ValueError(f'The given graph format ({graph_format}) is not correct. it should be "graphml" or gexf')

    return out_file


def compute_edge_metrics(context_graph: nx.Graph, gene_proportion_cutoff: float) -> None:
    """
    Compute various metrics on the edges of the context graph.

    :param context_graph: The context graph.
    :param gene_proportion_cutoff: The minimum proportion of shared genes between two features for their edge to be considered significant.
    """
    # compute jaccard on organism and on genes
    for f1, f2, data in context_graph.edges(data=True):
        data['jaccard_organism'] = len(data['organisms']) / len(set(f1.organisms) | set(f2.organisms))

        f1_gene_proportion = len(data['genes'][f1]) / f1.number_of_genes
        f2_gene_proportion = len(data['genes'][f2]) / f2.number_of_genes

        data['f1'] = f1.name
        data['f2'] = f2.name
        data['f1_jaccard_gene'] = f1_gene_proportion
        data['f2_jaccard_gene'] = f2_gene_proportion

        data[f'is_jaccard_gene_>_{gene_proportion_cutoff}'] = (f1_gene_proportion >= gene_proportion_cutoff) and (
                    f2_gene_proportion >= gene_proportion_cutoff)

        transitivity_counter = data['transitivity']

        mean_transitivity = sum(
            (transitivity * counter for transitivity, counter in transitivity_counter.items())) / sum(
            (counter for counter in transitivity_counter.values()))

        data['mean_transitivity'] = mean_transitivity

        # the following commented out lines are additional metrics that could be used

        # data['min_jaccard_organism'] = len(data['organisms'])/min(len(f1.organisms), len(f2.organisms))
        # data['max_jaccard_organism'] = len(data['organisms'])/max(len(f1.organisms), len(f2.organisms))
        # f1_gene_proportion_partial = len(data['genes'][f1])/len(context_graph.nodes[f1]['genes'])
        # f2_gene_proportion_partial = len(data['genes'][f2])/len(context_graph.nodes[f2]['genes'])
        # data[f'f1_jaccard_gene_partital'] = f1_gene_proportion_partial
        # data[f'f2_jaccard_gene_partital'] = f2_gene_proportion_partial


def add_edges_to_context_graph(context_graph: nx.Graph,
                               contig_genes: List[Gene],
                               contig_windows: List[Tuple[int, int]],
                               transitivity: int,
                               is_circular: bool):
    """
    Add edges to the context graph based on contig genes and windows.

    :param context_graph: The context graph to which edges will be added.
    :param contig_genes: An iterable of genes in the contig.
    :param contig_windows: A list of tuples representing the start and end positions of contig windows.
    :param transitivity: The number of next genes to consider when adding edges.
    :param is_circular: A boolean indicating if the contig is circular.

    """
    for window_start, window_end in contig_windows:
        for gene_index in range(window_start, window_end + 1):
            gene = contig_genes[gene_index]
            next_genes = get_n_next_genes_index(gene_index, next_genes_count=transitivity + 1,
                                                contig_size=len(contig_genes), is_circular=is_circular)
            next_genes = list(next_genes)

            for i, next_gene_index in enumerate(next_genes):
                # Check if the next gene is within the contig windows
                if not any(lower <= next_gene_index <= upper for (lower, upper) in contig_windows):
                    # next_gene_index is not in any range of genes in the context
                    # so it is ignored along with all following genes
                    break

                next_gene = contig_genes[next_gene_index]
                if next_gene.family == gene.family:
                    # If the next gene has the same family, the two genes refer to the same node
                    # so they are ignored
                    continue

                context_graph.add_edge(gene.family, next_gene.family)

                edge_dict = context_graph[gene.family][next_gene.family]

                if i == 0:
                    edge_dict['adjacent_family'] = True

                # Store information of the transitivity used to link the two genes:
                if "transitivity" not in edge_dict:
                    edge_dict['transitivity'] = {i: 0 for i in range(transitivity + 1)}
                edge_dict['transitivity'][i] += 1

                # Add node attributes
                node_gene_dict = context_graph.nodes[gene.family]
                next_gene_gene_dict = context_graph.nodes[next_gene.family]

                increment_attribute_counter(node_gene_dict, "genes_count")
                increment_attribute_counter(next_gene_gene_dict, "genes_count")

                add_val_to_dict_attribute(node_gene_dict, "genes", gene)
                add_val_to_dict_attribute(next_gene_gene_dict, "genes", next_gene)

                # Add edge attributes
                edge_dict = context_graph[gene.family][next_gene.family]
                try:
                    genes_edge_dict = edge_dict['genes']
                except Exception:
                    genes_edge_dict = {}
                    edge_dict['genes'] = genes_edge_dict

                add_val_to_dict_attribute(genes_edge_dict, gene.family, gene)
                add_val_to_dict_attribute(genes_edge_dict, next_gene.family, next_gene)

                add_val_to_dict_attribute(edge_dict, "organisms", gene.organism)

                increment_attribute_counter(edge_dict, "gene_pairs")

                assert gene.organism == next_gene.organism, (f"Gene of the same contig have a different organism. "
                                                             f"{gene.organism} and {next_gene.organism}")


def add_val_to_dict_attribute(attr_dict: dict, attribute_key, attribute_value):
    """
    Add an attribute value to a edge or node dictionary set.

    :param attr_dict: The dictionary containing the edge/node attributes.
    :param attribute_key: The key of the attribute.
    :param attribute_value: The value of the attribute to be added.

    """

    try:
        attr_dict[attribute_key].add(attribute_value)
    except KeyError:
        attr_dict[attribute_key] = {attribute_value}


def increment_attribute_counter(edge_dict: dict, key: Hashable):
    """
    Increment the counter for an edge/node attribute in the edge/node dictionary.

    :param edge_dict: The dictionary containing the attributes.
    :param key: The key of the attribute.

    """

    try:
        edge_dict[key] += 1
    except KeyError:
        edge_dict[key] = 1


def get_n_next_genes_index(current_index: int, next_genes_count: int,
                           contig_size: int, is_circular: bool = False) -> Iterator[int]:
    """
    Generate the indices of the next genes based on the current index and contig properties.

    :param current_index: The index of the current gene.
    :param next_genes_count: The number of next genes to consider.
    :param contig_size: The total number of genes in the contig.
    :param is_circular: Flag indicating whether the contig is circular (default: False).
    :return: An iterator yielding the indices of the next genes.

    Raises:
    - IndexError: If the current index is out of range for the given contig size.

    """

    # Check if the current index is out of range
    if current_index >= contig_size:
        raise IndexError(f'current gene index is out of range. '
                         f"Contig has {contig_size} genes while the given gene index is {current_index}")
    if is_circular:
        next_genes = chain(range(current_index + 1, contig_size), range(0, current_index))
    else:
        next_genes = range(current_index + 1, contig_size)

    for i, next_gene_index in enumerate(next_genes):
        if i == next_genes_count:
            break
        yield next_gene_index


def get_contig_to_genes(gene_families: Iterable[GeneFamily]) -> Dict[Contig, Set[Gene]]:
    """
    Group genes from specified gene families by contig.

    :param gene_families: An iterable of gene families object.
    
    :return: A dictionary mapping contigs to sets of genes.
    """

    contig_to_genes_of_interest = defaultdict(set)
    for gene_family in gene_families:
        for gene in gene_family.genes:
            contig = gene.contig
            contig_to_genes_of_interest[contig].add(gene)
    return contig_to_genes_of_interest


def compute_gene_context_graph(families: Iterable[GeneFamily], transitive: int = 4, window_size: int = 0,
                               disable_bar: bool = False) -> nx.Graph:
    """
    Construct the graph of gene contexts between families of the pangenome.

    :param families: An iterable of gene families.
    :param transitive: Size of the transitive closure used to build the graph.
    :param window_size: Size of the window for extracting gene contexts (default: 0).
    :param disable_bar: Flag to disable the progress bar (default: False).

    :return: The constructed gene context graph.
    """

    context_graph = nx.Graph()

    contig_to_genes_of_interest = get_contig_to_genes(families)

    for contig, genes_of_interest in tqdm(contig_to_genes_of_interest.items(), unit="contig",
                                          total=len(contig_to_genes_of_interest), disable=disable_bar):
        genes_count = contig.number_of_genes

        genes_of_interest_positions = [g.position for g in genes_of_interest]

        contig_windows = extract_contig_window(genes_count, genes_of_interest_positions,
                                               window_size=window_size, is_circular=contig.is_circular)

        add_edges_to_context_graph(context_graph,
                                   contig.get_genes(),
                                   contig_windows,
                                   transitive,
                                   contig.is_circular)
    return context_graph


def fam_to_seq(seq_to_pan: dict) -> dict:
    """
    Create a dictionary with gene families as keys and list of sequences id as values

    :param seq_to_pan: Dictionary storing the sequence ids as keys and the gene families
                       to which they are assigned as values

    :return: Dictionary reversed
    """

    fam_2_seq = {}
    for sequence, family in seq_to_pan.items():
        if family.ID in fam_2_seq.keys():
            fam_2_seq[family.ID].append(sequence)
        else:
            fam_2_seq[family.ID] = [sequence]
    return fam_2_seq


def export_context_to_dataframe(gene_contexts: set, fam2seq: Dict[str, int],
                                families_of_interest: Set[GeneFamily], output: Path):
    """
    Export the results into dataFrame

    :param gene_contexts: connected components found in the pangenome
    :param fam2seq: Dictionary with gene families ID as keys and list of sequence ids as values
    :param families_of_interest: families of interest that are at the origine of the context.
    :param output: output path
    """

    lines = []
    for gene_context in gene_contexts:
        for family in gene_context.families:
            if fam2seq.get(family) is None:
                sequence_id = None
            else:
                sequence_id = ','.join(fam2seq.get(family))

            family_info = {"GeneContext ID": gene_context.ID,
                           "Gene family name": family.name,
                           "Sequence ID": sequence_id,
                           "Nb Genomes": family.number_of_organisms,
                           "Partition": family.named_partition,
                           "Target family": family in families_of_interest}
            
            lines.append(family_info)

    df = pd.DataFrame(lines).set_index("GeneContext ID")

    df = df.sort_values(["GeneContext ID", "Sequence ID"], na_position='last')

    df.to_csv(output, sep="\t", na_rep='NA')

    logging.getLogger().debug(f"detected gene context(s) are listed in: '{output}'")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """

    if not any([args.sequences, args.family]):
        raise Exception("At least one of --sequences or --family option must be given")

    mk_outdir(args.output, args.force)

    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    # check statuses and load info
    if args.sequences is not None and pangenome.status["geneFamilySequences"] not in ["inFile", "Loaded", "Computed"]:
        raise Exception("Cannot use this function as your pangenome does not have gene families representatives "
                        "associated to it. For now this works only if the clustering has been made by PPanGGOLiN.")

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=args.disable_prog_bar)

    search_gene_context_in_pangenome(pangenome=pangenome, output=args.output, tmpdir=args.tmpdir,
                                     sequence_file=args.sequences, families=args.family, transitive=args.transitive,
                                     identity=args.identity, coverage=args.coverage, use_representatives=args.fast,
                                     jaccard_threshold=args.jaccard,
                                     window_size=args.window_size,
                                     no_defrag=args.no_defrag, cpu=args.cpu, disable_bar=args.disable_prog_bar,
                                     graph_format=args.graph_format,
                                     translation_table=args.translation_table, keep_tmp=args.keep_tmp)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """

    parser = sub_parser.add_parser("context", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_context(parser)
    return parser


def parser_context(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of context command

    :param parser: parser for align argument
    """

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome.h5 file")
    required.add_argument('-o', '--output', required=False, type=Path,
                          default="ppanggolin_context" + time.strftime("_DATE%Y-%m-%d_HOUR%H.%M.%S",
                                                                       time.localtime()) + "_PID" + str(os.getpid()),
                          help="Output directory where the file(s) will be written")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('-S', '--sequences', required=False, type=Path,
                        help="Fasta file with the sequences of interest")
    onereq.add_argument('-F', '--family', required=False, type=Path,
                        help="List of family IDs of interest from the pangenome")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    optional.add_argument("--fast", required=False, action="store_true",
                          help="Use representative sequences of gene families for input gene alignment. "
                               "This option is recommended for faster processing but may be less sensitive. "
                               "By default, all pangenome genes are used for alignment. "
                               "This argument makes sense only when --sequence is provided.")
    optional.add_argument('--identity', required=False, type=float, default=0.8,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("--translation_table", required=False, default="11",
                          help="The translation table (genetic code) to use when the input sequences are nucleotide "
                               "sequences. ")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. This indicates the number of "
                               "non related genes allowed in-between two related genes. Increasing it will improve "
                               "precision but lower sensitivity a little.")
    optional.add_argument("-w", "--window_size", required=False, type=int, default=5,
                          help="Number of neighboring genes that are considered on each side of "
                               "a gene of interest when searching for conserved genomic contexts.")

    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
                               "will improve precision but lower sensitivity a lot.")
    optional.add_argument('--graph_format', help="Format of the context graph. Can be gexf or graphml.",
                          default='graphml', choices=['gexf', 'graphml'])
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int,
                          help="Number of available cpus")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_context(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
