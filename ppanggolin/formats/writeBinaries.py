#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from collections import Counter, defaultdict
import statistics
from typing import Tuple
import pkg_resources

# installed libraries
from tqdm import tqdm
import tables
from gmpy2 import popcount

#local libraries
from ppanggolin.pangenome import Pangenome


def gene_desc(org_len, contig_len, id_len, type_len, name_len, product_len, max_local_id) -> dict:
    """
    Create a table to save gene description

    :param org_len: Maximum size of organism
    :param contig_len: Maximum size of contigs
    :param id_len: Maximum size of gene ID
    :param type_len: Maximum size of gene Type
    :param name_len: Maximum size of gene name
    :param product_len: Maximum size of gene product
    :param max_local_id: Maximum size of gene local identifier

    :return: Formatted table
    """
    return {
        'organism': tables.StringCol(itemsize=org_len),
        "contig": {
            'name': tables.StringCol(itemsize=contig_len),
            "is_circular": tables.BoolCol(dflt=False)
        },
        "gene": {
            'ID': tables.StringCol(itemsize=id_len),
            'start': tables.UInt64Col(),
            'stop': tables.UInt64Col(),
            'strand': tables.StringCol(itemsize=1),
            'type': tables.StringCol(itemsize=type_len),
            'position': tables.UInt32Col(),
            'name': tables.StringCol(itemsize=name_len),
            'product': tables.StringCol(itemsize=product_len),
            'genetic_code': tables.UInt32Col(dflt=11),
            'is_fragment': tables.BoolCol(dflt=False),
            'local': tables.StringCol(itemsize=max_local_id)
        }
    }


def get_max_len_annotations(pangenome: Pangenome) -> Tuple[int, int, int, int, int, int, int]:
    """
    Get the maximum size of each annotation information to optimize saving

    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_org_len = 1
    max_contig_len = 1
    max_gene_id_len = 1
    max_name_len = 1
    max_product_len = 1
    max_type_len = 1
    max_local_id = 1
    for org in pangenome.organisms:
        if len(org.name) > max_org_len:
            max_org_len = len(org.name)
        for contig in org.contigs:
            if len(contig.name) > max_contig_len:
                max_contig_len = len(contig.name)
            for gene in contig.genes:
                if len(gene.ID) > max_gene_id_len:
                    max_gene_id_len = len(gene.ID)
                if len(gene.name) > max_name_len:
                    max_name_len = len(gene.name)
                if len(gene.product) > max_product_len:
                    max_product_len = len(gene.product)
                if len(gene.type) > max_type_len:
                    max_type_len = len(gene.type)
                if len(gene.local_identifier) > max_local_id:
                    max_local_id = len(gene.local_identifier)
            for gene in contig.RNAs:
                if len(gene.ID) > max_gene_id_len:
                    max_gene_id_len = len(gene.ID)
                if len(gene.name) > max_name_len:
                    max_name_len = len(gene.name)
                if len(gene.product) > max_product_len:
                    max_product_len = len(gene.product)
                if len(gene.type) > max_type_len:
                    max_type_len = len(gene.type)
                if len(gene.local_identifier) > max_local_id:
                    max_local_id = len(gene.local_identifier)

    return max_org_len, max_contig_len, max_gene_id_len, max_type_len, max_name_len, max_product_len, max_local_id

class geneTableElement:
    def __init__(self, start, stop, strand, genetype, position, name, product, is_fragment, genetic_code):
        """stores potentially non unique descriptive info of a gene"""
        self.start = start
        self.stop = stop
        self.strand = strand
        self.type=genetype
        self.position=position
        self.name=name
        self.product=product
        self.is_fragment=is_fragment
        self.genetic_code = genetic_code

    def get_element(self):
        return (self.start, self.stop, self.strand, self.type, self.position, self.name, self.product, self.is_fragment, self.genetic_code)

def write_annotations(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Function writing all the pangenome annotations

    :param pangenome: Annotated pangenome
    :param h5f: Pangenome HDF5 file
    :param disable_bar: Alow to disable progress bar
    """
    annotation = h5f.create_group("/", "annotations", "Annotations of the pangenome organisms")
    gene_table = h5f.create_table(annotation, "genes", gene_desc(*get_max_len_annotations(pangenome)),
                                  expectedrows=len(pangenome.genes))
    geneTabs = set()

    gene_row = gene_table.row
    for org in tqdm(pangenome.organisms, total=pangenome.number_of_organisms(), unit="genome", disable=disable_bar):
        for contig in org.contigs:
            for gene in contig.genes:
                gene_row["organism"] = org.name
                gene_row["contig/name"] = contig.name
                gene_row["contig/is_circular"] = contig.is_circular  # this should be somewhere else.
                gene_row["gene/ID"] = gene.ID
                gene_row["gene/start"] = gene.start
                gene_row["gene/stop"] = gene.stop
                gene_row["gene/strand"] = gene.strand
                gene_row["gene/type"] = gene.type
                gene_row["gene/position"] = gene.position
                gene_row["gene/name"] = gene.name
                gene_row["gene/product"] = gene.product
                gene_row["gene/is_fragment"] = gene.is_fragment
                gene_row["gene/genetic_code"] = gene.genetic_code
                gene_row["gene/local"] = gene.local_identifier
                gene_row.append()
                geneTabs.add((gene.start, gene.stop, gene.strand, gene.type, gene.position, gene.name, gene.product, gene.is_fragment, gene.genetic_code))
            for rna in contig.RNAs:
                gene_row["organism"] = org.name
                gene_row["contig/name"] = contig.name
                gene_row["contig/is_circular"] = contig.is_circular  # this should be somewhere else.
                gene_row["gene/ID"] = rna.ID
                gene_row["gene/start"] = rna.start
                gene_row["gene/stop"] = rna.stop
                gene_row["gene/strand"] = rna.strand
                gene_row["gene/type"] = rna.type
                gene_row["gene/name"] = rna.name
                gene_row["gene/product"] = rna.product
                gene_row["gene/is_fragment"] = rna.is_fragment
                gene_row.append()
    gene_table.flush()

    logging.getLogger().info(f"would get {len(geneTabs)} elements instead of {len(pangenome.genes)}")

    #for val in sorted(geneTabs, key=lambda x:x[0]):
    #    print(val)


def get_gene_sequences_len(pangenome: Pangenome) -> Tuple[int, int]:
    """
    Get the maximum size of gene sequences to optimize disk space

    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_gene_id_len = 1
    max_gene_type = 1
    for gene in pangenome.genes:
        if len(gene.ID) > max_gene_id_len:
            max_gene_id_len = len(gene.ID)
        if len(gene.type) > max_gene_type:
            max_gene_type = len(gene.type)
    return max_gene_id_len, max_gene_type


def gene_sequences_desc(gene_id_len, gene_type_len) -> dict:
    """
    Create table to save gene sequences

    :param gene_id_len: Maximum size of gene sequence identifier
    :param gene_type_len: Maximum size of gene type

    :return: Formated table
    """
    return {
        "gene": tables.StringCol(itemsize=gene_id_len),
        "seqid": tables.UInt32Col(),
        "type": tables.StringCol(itemsize=gene_type_len)
    }

def get_sequence_len(pangenome: Pangenome) -> int:
    """
    Get the maximum size of gene sequences to optimize disk space

    :param pangenome: Annotated pangenome
    :return: maximum size of each annotation
    """
    max_seq_len = 1
    for gene in pangenome.genes:
        if len(gene.dna) > max_seq_len:
            max_seq_len = len(gene.dna)
    return max_seq_len

def sequence_desc(max_seq_len: int) -> dict:
    """
    Table description to save sequences
    :param max_seq_len: Maximum size of gene type

    :return: Formated table
    """
    return {
        "seqid": tables.UInt32Col(),
        "dna": tables.StringCol(itemsize=max_seq_len)
    }

def write_gene_sequences(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Function writing all the pangenome gene sequences

    :param pangenome: Pangenome with gene sequences
    :param h5f: Pangenome HDF5 file without sequences
    :param disable_bar: Disable progress bar
    """
    gene_seq = h5f.create_table("/", "geneSequences", gene_sequences_desc(*get_gene_sequences_len(pangenome)),
                                expectedrows=len(pangenome.genes))
    #process sequences to save them only once
    seq2seqid = {}
    id_counter = 0
    gene_row = gene_seq.row
    for gene in tqdm(pangenome.genes, total=pangenome.number_of_gene(), unit="gene", disable=disable_bar):
        curr_seq_id = seq2seqid.get(gene.dna)
        if curr_seq_id is None:
            curr_seq_id = id_counter
            seq2seqid[gene.dna] = id_counter
            id_counter+=1
        gene_row["gene"] = gene.ID
        gene_row["seqid"] = curr_seq_id
        gene_row["type"] = gene.type
        gene_row.append()
    gene_seq.flush()

    seq_table = h5f.create_table("/","sequences", sequence_desc(get_sequence_len(pangenome)),
                                 expectedrows=len(seq2seqid))

    seq_row = seq_table.row
    for seq, seqid in seq2seqid.items():
        seq_row["dna"] = seq
        seq_row["seqid"] = seqid
        seq_row.append()
    seq_table.flush()


def gene_fam_desc(max_name_len: int, max_sequence_length: int, max_part_len: int) -> dict:
    """
    Create a formated table for gene families description

    :param max_name_len: Maximum size of gene family name
    :param max_sequence_length: Maximum size of gene family representing gene sequences
    :param max_part_len: Maximum size of gene family partition

    :return: Formated table
    """
    return {
        "name": tables.StringCol(itemsize=max_name_len),
        "protein": tables.StringCol(itemsize=max_sequence_length),
        "partition": tables.StringCol(itemsize=max_part_len)
    }


def get_gene_fam_len(pangenome: Pangenome) -> Tuple[int, int, int]:
    """
    Get maximum size of gene families information

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_fam_name_len = 1
    max_gene_fam_seq_len = 1
    max_part_len = 1
    for genefam in pangenome.gene_families:
        if len(genefam.sequence) > max_gene_fam_seq_len:
            max_gene_fam_seq_len = len(genefam.sequence)
        if len(genefam.name) > max_gene_fam_name_len:
            max_gene_fam_name_len = len(genefam.name)
        if len(genefam.partition) > max_part_len:
            max_part_len = len(genefam.partition)
    return max_gene_fam_name_len, max_gene_fam_seq_len, max_part_len


def write_gene_fam_info(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family

    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param force: force to write information if precedent information exist
    :param disable_bar: Disable progress bar
    """
    if '/geneFamiliesInfo' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')  # erasing the table, and rewriting a new one.
    gene_fam_seq = h5f.create_table("/", "geneFamiliesInfo", gene_fam_desc(*get_gene_fam_len(pangenome)),
                                    expectedrows=len(pangenome.gene_families))

    row = gene_fam_seq.row
    for fam in tqdm(pangenome.gene_families, total=pangenome.number_of_gene_families(),
                    unit="gene family", disable=disable_bar):
        row["name"] = fam.name
        row["protein"] = fam.sequence
        row["partition"] = fam.partition
        row.append()
    gene_fam_seq.flush()


def gene_to_fam_desc(gene_fam_name_len: int, gene_id_len: int) -> dict:
    """
    Create a formated table for gene in gene families information

    :param gene_fam_name_len: Maximum size of gene family names
    :param gene_id_len: Maximum sizez of gene identifier

    :return: formated table
    """
    return {
        "geneFam": tables.StringCol(itemsize=gene_fam_name_len),
        "gene": tables.StringCol(itemsize=gene_id_len)
    }


def get_gene_to_fam_len(pangenome: Pangenome):
    """
    Get maximum size of gene in gene families information

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_fam_name = 1
    max_gene_id = 1
    for geneFam in pangenome.gene_families:
        if len(geneFam.name) > max_gene_fam_name:
            max_gene_fam_name = len(geneFam.name)
        for gene in geneFam.genes:
            if len(gene.ID) > max_gene_id:
                max_gene_id = len(gene.ID)
    return max_gene_fam_name, max_gene_id


def write_gene_families(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Function writing all the pangenome gene families

    :param pangenome: pangenome with gene families computed
    :param h5f: HDF5 file to save pangenome with gene families
    :param force: Force to write gene families in hdf5 file if there is already gene families
    :param disable_bar: Disable progress bar
    """
    if '/gene_families' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'gene_families')  # erasing the table, and rewriting a new one.
    gene_families = h5f.create_table("/", "geneFamilies", gene_to_fam_desc(*get_gene_to_fam_len(pangenome)))
    gene_row = gene_families.row
    for geneFam in tqdm(pangenome.gene_families, total=pangenome.number_of_gene_families(), unit="gene family",
                        disable=disable_bar):
        for gene in geneFam.genes:
            gene_row["gene"] = gene.ID
            gene_row["geneFam"] = geneFam.name
            gene_row.append()
    gene_families.flush()


def graph_desc(max_gene_id_len):
    """
    Create a formated table for pangenome graph

    :param max_gene_id_len: Maximum size of gene id

    :return: formated table
    """
    return {
        'geneTarget': tables.StringCol(itemsize=max_gene_id_len),
        'geneSource': tables.StringCol(itemsize=max_gene_id_len)
    }


def get_gene_id_len(pangenome: Pangenome) -> int:
    """
    Get maximum size of gene id in pangenome graph

    :param pangenome: Pangenome with graph computed

    :return: Maximum size of gene id
    """
    max_gene_len = 1
    for gene in pangenome.genes:
        if len(gene.ID) > max_gene_len:
            max_gene_len = len(gene.ID)
    return max_gene_len


def write_graph(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Function writing the pangenome graph

    :param pangenome: pangenome with graph computed
    :param h5f: HDF5 file to save pangenome graph
    :param force: Force to write graph in hdf5 file if there is already one
    :param disable_bar: Disable progress bar
    """
    # TODO if we want to be able to read the graph without reading the annotations (because it's one of the most time
    #  consumming parts to read), it might be good to add the organism name in the table here.
    #  for now, forcing the read of annotations.
    if '/edges' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/", "edges")
    edge_table = h5f.create_table("/", "edges", graph_desc(get_gene_id_len(pangenome)),
                                  expectedrows=len(pangenome.edges))
    edge_row = edge_table.row
    for edge in tqdm(pangenome.edges, total=pangenome.number_of_edge(), unit="edge", disable=disable_bar):
        for genePairs in edge.organisms.values():
            for gene1, gene2 in genePairs:
                edge_row["geneTarget"] = gene1.ID
                edge_row["geneSource"] = gene2.ID
                edge_row.append()
    edge_table.flush()


def rgp_desc(max_rgp_len, max_gene_len):
    """
    Create a formated table for region of genomic plasticity

    :param max_rgp_len: Maximum size of RGP
    :param max_gene_len: Maximum sizez of gene

    :return: formated table
    """
    return {
        'RGP': tables.StringCol(itemsize=max_rgp_len),
        'gene': tables.StringCol(itemsize=max_gene_len)
    }


def get_rgp_len(pangenome: Pangenome) -> Tuple[int, int]:
    """
    Get maximum size of region of genomic plasticity and gene

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_gene_len = 1
    max_rgp_len = 1
    for region in pangenome.regions:
        for gene in region.genes:
            if len(gene.ID) > max_gene_len:
                max_gene_len = len(gene.ID)
        if len(region.name) > max_rgp_len:
            max_rgp_len = len(region.name)
    return max_rgp_len, max_gene_len


def write_rgp(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Function writing all the region of genomic plasticity in pangenome

    :param pangenome: pangenome with RGP computed
    :param h5f: HDF5 file to save pangenome with RGP
    :param force: Force to write gene families in hdf5 file if there is already RGP
    :param disable_bar: Disable progress bar
    """
    if '/RGP' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computer RGP")
        h5f.remove_node('/', 'RGP')

    rgp_table = h5f.create_table('/', 'RGP', rgp_desc(*get_rgp_len(pangenome)),
                                 expectedrows=sum([len(region.genes) for region in pangenome.regions]))
    rgp_row = rgp_table.row
    for region in tqdm(pangenome.regions, total=pangenome.number_of_rgp(), unit="region", disable=disable_bar):
        for gene in region.genes:
            rgp_row["RGP"] = region.name
            rgp_row["gene"] = gene.ID
            rgp_row.append()
    rgp_table.flush()


def spot_desc(max_rgp_len):
    """
    Create a formated table for hotspot

    :param max_rgp_len: Maximum size of RGP

    :return: formated table
    """
    return {
        'spot': tables.UInt32Col(),
        'RGP': tables.StringCol(itemsize=max_rgp_len)
    }


def get_spot_desc(pangenome: Pangenome) -> int:
    """
    Get maximum size of region of genomic plasticity in hotspot

    :param pangenome: Pangenome with gene families computed

    :return: Maximum size of each element
    """
    max_rgp_len = 1
    for spot in pangenome.spots:
        for region in spot.regions:
            if len(region.name) > max_rgp_len:
                max_rgp_len = len(region.name)
    return max_rgp_len


def write_spots(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Function writing all the pangenome hotspot

    :param pangenome: pangenome with spot computed
    :param h5f: HDF5 file to save pangenome with spot
    :param force: Force to write gene families in hdf5 file if there is already spot
    :param disable_bar: Disable progress bar
    """
    if '/spots' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed spots")
        h5f.remove_node("/", "spots")

    spot_table = h5f.create_table("/", "spots", spot_desc(get_spot_desc(pangenome)),
                                  expectedrows=sum([len(spot.regions) for spot in pangenome.spots]))
    spot_row = spot_table.row
    for spot in tqdm(pangenome.spots, total=pangenome.number_of_spots(), unit="spot", disable=disable_bar):
        for region in spot.regions:
            spot_row["spot"] = spot.ID
            spot_row["RGP"] = region.name
            spot_row.append()
    spot_table.flush()


def mod_desc(gene_fam_name_len):
    """
    Create a formated table for hotspot

    :param gene_fam_name_len: Maximum size of gene families name

    :return: formated table
    """
    return {
        "geneFam": tables.StringCol(itemsize=gene_fam_name_len),
        "module": tables.UInt32Col(),
    }


def get_mod_desc(pangenome: Pangenome) -> int:
    """
    Get maximum size of gene families name in modules

    :param pangenome: Pangenome with modules computed

    :return: Maximum size of each element
    """
    max_fam_len = 1
    for mod in pangenome.modules:
        for fam in mod.families:
            if len(fam.name) > max_fam_len:
                max_fam_len = len(fam.name)
    return max_fam_len


def write_modules(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Function writing all the pangenome modules

    :param pangenome: pangenome with spot computed
    :param h5f: HDF5 file to save pangenome with spot
    :param force: Force to write gene families in hdf5 file if there is already spot
    :param disable_bar: Disable progress bar
    """
    if '/modules' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed modules")
        h5f.remove_node("/", "modules")

    mod_table = h5f.create_table('/', 'modules', mod_desc(get_mod_desc(pangenome)),
                                 expectedrows=sum([len(mod.families) for mod in pangenome.modules]))
    mod_row = mod_table.row

    for mod in tqdm(pangenome.modules, total=pangenome.number_of_modules(), unit="modules", disable=disable_bar):
        for fam in mod.families:
            mod_row["geneFam"] = fam.name
            mod_row["module"] = mod.ID
            mod_row.append()
    mod_table.flush()


def write_status(pangenome: Pangenome, h5f: tables.File):
    """
    Write pangenome status in HDF5 file

    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    """
    if "/status" in h5f:  # if statuses are already written
        status_group = h5f.root.status
    else:  # else create the status group.
        status_group = h5f.create_group("/", "status", "Statuses of the pangenome content")
    status_group._v_attrs.genomesAnnotated = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded",
                                                                                              "inFile"] else False
    status_group._v_attrs.geneSequences = True if pangenome.status["geneSequences"] in ["Computed", "Loaded",
                                                                                        "inFile"] else False
    status_group._v_attrs.genesClustered = True if pangenome.status["genesClustered"] in ["Computed", "Loaded",
                                                                                          "inFile"] else False
    status_group._v_attrs.geneFamilySequences = True if pangenome.status["geneFamilySequences"] in ["Computed",
                                                                                                    "Loaded",
                                                                                                    "inFile"] else False
    status_group._v_attrs.NeighborsGraph = True if pangenome.status["neighborsGraph"] in ["Computed", "Loaded",
                                                                                          "inFile"] else False
    status_group._v_attrs.Partitioned = True if pangenome.status["partitioned"] in ["Computed", "Loaded",
                                                                                    "inFile"] else False
    status_group._v_attrs.defragmented = True if pangenome.status["defragmented"] in ["Computed", "Loaded",
                                                                                      "inFile"] else False
    status_group._v_attrs.predictedRGP = True if pangenome.status["predictedRGP"] in ["Computed", "Loaded",
                                                                                      "inFile"] else False
    status_group._v_attrs.spots = True if pangenome.status["spots"] in ["Computed", "Loaded", "inFile"] else False
    status_group._v_attrs.modules = True if pangenome.status["modules"] in ["Computed", "Loaded", "inFile"] else False
    status_group._v_attrs.version = pkg_resources.get_distribution("ppanggolin").version


def write_info(pangenome: Pangenome, h5f: tables.File):
    """
    Writes information and numbers to be eventually called with the 'info' submodule

    :param pangenome: Pangenome object with some information computed
    :param h5f: Pangenome file to save information
    """

    def getmean(arg: iter) -> float:
        """ Compute the mean of arguments if exist 0 else

        :param arg: list of values

        :return: return the mean
        """
        return 0 if len(arg) == 0 else round(statistics.mean(arg), 2)

    def getstdev(arg: iter) -> float:
        """ Compute the standard deviation of arguments if exist 0 else

        :param arg: list of values

        :return: return the sd
        """
        return 0 if len(arg) <= 1 else round(statistics.stdev(arg), 2)

    def getmax(arg: iter) -> float:
        """ Get the maximum of arguments if exist 0 else

        :param arg: list of values

        :return: return the maximum
        """
        return 0 if len(arg) == 0 else round(max(arg), 2)

    def getmin(arg: iter) -> float:
        """ Get the minimum of arguments if exist 0 else

        :param arg: list of values

        :return: return the minimum
        """
        return 0 if len(arg) == 0 else round(min(arg), 2)

    if "/info" in h5f:
        info_group = h5f.root.info
    else:
        info_group = h5f.create_group("/", "info", "Informations about the pangenome content")
    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfGenes = len(pangenome.genes)
        info_group._v_attrs.numberOfOrganisms = len(pangenome.organisms)
    if pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfClusters = len(pangenome.gene_families)
    if pangenome.status["neighborsGraph"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfEdges = len(pangenome.edges)
    if pangenome.status["partitioned"] in ["Computed", "Loaded"]:
        named_part_counter = Counter()
        subpart_counter = Counter()
        part_distribs = defaultdict(list)
        part_set = set()
        for fam in pangenome.gene_families:
            named_part_counter[fam.named_partition] += 1
            part_distribs[fam.named_partition].append(len(fam.organisms) / len(pangenome.organisms))
            if fam.named_partition == "shell":
                subpart_counter[fam.partition] += 1
            if fam.partition != "S_":
                part_set.add(fam.partition)

        info_group._v_attrs.numberOfPersistent = named_part_counter["persistent"]
        info_group._v_attrs.persistentStats = {"min": getmin(part_distribs["persistent"]),
                                               "max": getmax(part_distribs["persistent"]),
                                               "sd": getstdev(part_distribs["persistent"]),
                                               "mean": getmean(part_distribs["persistent"])}
        info_group._v_attrs.numberOfShell = named_part_counter["shell"]
        info_group._v_attrs.shellStats = {"min": getmin(part_distribs["shell"]), "max": getmax(part_distribs["shell"]),
                                          "sd": getstdev(part_distribs["shell"]),
                                          "mean": getmean(part_distribs["shell"])}
        info_group._v_attrs.numberOfCloud = named_part_counter["cloud"]
        info_group._v_attrs.cloudStats = {"min": getmin(part_distribs["cloud"]), "max": getmax(part_distribs["cloud"]),
                                          "sd": getstdev(part_distribs["cloud"]),
                                          "mean": getmean(part_distribs["cloud"])}
        info_group._v_attrs.numberOfPartitions = len(part_set)
        info_group._v_attrs.numberOfSubpartitions = subpart_counter
    if pangenome.status["predictedRGP"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfRGP = len(pangenome.regions)
    if pangenome.status["spots"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfSpots = len(pangenome.spots)
    if pangenome.status["modules"] in ["Computed", "Loaded"]:
        info_group._v_attrs.numberOfModules = len(pangenome.modules)
        info_group._v_attrs.numberOfFamiliesInModules = sum([len(mod.families) for mod in pangenome.modules])

    info_group._v_attrs.parameters = pangenome.parameters  # saving the pangenome parameters


def write_info_modules(pangenome: Pangenome, h5f: tables.File):
    """
    Writes more information about modules if computed by metrics subpackage

    :param pangenome: Pangenome object with some information computed
    :param h5f: Pangenome file to save information
    """

    def getmean(arg: iter) -> float:
        """ Compute the mean of arguments if exist 0 else

        :param arg: list of values

        :return: return the mean
        """
        return 0 if len(arg) == 0 else round(statistics.mean(arg), 2)

    def getstdev(arg: iter) -> float:
        """ Compute the standard deviation of arguments if exist 0 else

        :param arg: list of values

        :return: return the sd
        """
        return 0 if len(arg) <= 1 else round(statistics.stdev(arg), 2)

    def getmax(arg: iter) -> float:
        """ Get the maximum of arguments if exist 0 else

        :param arg: list of values

        :return: return the maximum
        """
        return 0 if len(arg) == 0 else round(max(arg), 2)

    def getmin(arg: iter) -> float:
        """ Get the minimum of arguments if exist 0 else

        :param arg: list of values

        :return: return the minimum
        """
        return 0 if len(arg) == 0 else round(min(arg), 2)

    if "/info" not in h5f:
        write_info(pangenome, h5f)
    info_group = h5f.root.info

    if pangenome.status["modules"] in ["Computed", "Loaded"]:
        def part_spec(part: str) -> list:
            """
            Get the list of module for a specific partition of pangenome

            :param part: pangenome partition name

            :return: list of module specific to partition
            """
            pangenome.compute_mod_bitarrays(part)
            return [popcount(module.bitarray) for module in pangenome.modules]

        mod_fam = [len(module.families) for module in pangenome.modules]
        info_group._v_attrs.StatOfFamiliesInModules = {"min": getmin(mod_fam),
                                                       "max": getmax(mod_fam),
                                                       "sd": getstdev(mod_fam),
                                                       "mean": getmean(mod_fam)}
        spec_pers = part_spec(part='persistent')
        spec_shell = part_spec(part='shell')
        spec_cloud = part_spec(part='cloud')
        info_group._v_attrs.PersistentSpecInModules = {"percent": round((sum(spec_pers) / sum(mod_fam)) * 100, 2),
                                                       "min": getmin(spec_pers),
                                                       "max": getmax(spec_pers),
                                                       "sd": getstdev(spec_pers),
                                                       "mean": getmean(spec_pers)}
        info_group._v_attrs.ShellSpecInModules = {"percent": round((sum(spec_shell) / sum(mod_fam)) * 100, 2),
                                                  "min": getmin(spec_shell),
                                                  "max": getmax(spec_shell),
                                                  "sd": getstdev(spec_shell),
                                                  "mean": getmean(spec_shell)}
        info_group._v_attrs.CloudSpecInModules = {"percent": round((sum(spec_cloud) / sum(mod_fam)) * 100, 2),
                                                  "min": getmin(spec_cloud),
                                                  "max": getmax(spec_cloud),
                                                  "sd": getstdev(spec_cloud),
                                                  "mean": getmean(spec_cloud)}
    else:
        raise Exception("Modules were not computed in your pangenome. Please see the module subcommand.")


def update_gene_fam_partition(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Update the gene families table with partition information

    :param pangenome: Partitioned pangenome
    :param h5f: HDF5 file with gene families
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger().info("Updating gene families with partition information")
    table = h5f.root.geneFamiliesInfo
    for row in tqdm(table, total=table.nrows, unit="gene family", disable=disable_bar):
        row["partition"] = pangenome.get_gene_family(row["name"].decode()).partition
        row.update()


def update_gene_fragments(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Updates the annotation table with the fragmentation information from the defrag pipeline

    :param pangenome: Annotated pangenome
    :param h5f: HDF5 pangenome file
    :param disable_bar: Allow to disable progress bar
    """
    logging.getLogger().info("Updating annotations with fragment information")
    table = h5f.root.annotations.genes
    for row in tqdm(table, total=table.nrows, unit="gene", disable=disable_bar):
        if row['gene/type'].decode() == 'CDS':
            row['gene/is_fragment'] = pangenome.get_gene(row['gene/ID'].decode()).is_fragment
            row.update()
    table.flush()


def erase_pangenome(pangenome: Pangenome, graph: bool = False, gene_families: bool = False, partition: bool = False,
                    rgp: bool = False, spots: bool = False, modules: bool = False):
    """
    Erases tables from a pangenome .h5 file

    :param pangenome: Pangenome
    :param graph: remove graph information
    :param gene_families: remove gene families information
    :param partition: remove partition information
    :param rgp: remove rgp information
    :param spots: remove spots information
    :param modules: remove modules information
    """

    h5f = tables.open_file(pangenome.file, "a")
    status_group = h5f.root.status
    info_group = h5f.root.info

    if '/edges' in h5f and (graph or gene_families):
        logging.getLogger().info("Erasing the formerly computed edges")
        h5f.remove_node("/", "edges")
        status_group._v_attrs.NeighborsGraph = False
        pangenome.status["neighborsGraph"] = "No"
        h5f.del_node_attr(info_group, "numberOfEdges")
    if '/gene_families' in h5f and gene_families:
        logging.getLogger().info("Erasing the formerly computed gene family to gene associations...")
        h5f.remove_node('/', 'gene_families')  # erasing the table, and rewriting a new one.
        pangenome.status["defragmented"] = "No"
        pangenome.status["genesClustered"] = "No"
        status_group._v_attrs.defragmented = False
        status_group._v_attrs.genesClustered = False

        h5f.del_node_attr(info_group, "numberOfClusters")

    if '/geneFamiliesInfo' in h5f and gene_families:
        logging.getLogger().info("Erasing the formerly computed gene family representative sequences...")
        h5f.remove_node('/', 'geneFamiliesInfo')  # erasing the table, and rewriting a new one.
        pangenome.status["geneFamilySequences"] = "No"
        status_group._v_attrs.geneFamilySequences = False
        if partition:
            logging.getLogger().info("Erasing former partitions...")
            pangenome.status["partitioned"] = "No"
            status_group._v_attrs.Partitioned = False
            if 'Partitionned' in status_group._v_attrs._f_list():
                status_group._v_attrs.Partitionned = False

            h5f.del_node_attr(info_group, "numberOfPersistent")
            h5f.del_node_attr(info_group, "persistentStats")
            h5f.del_node_attr(info_group, "numberOfShell")
            h5f.del_node_attr(info_group, "shellStats")
            h5f.del_node_attr(info_group, "numberOfCloud")
            h5f.del_node_attr(info_group, "cloudStats")
            h5f.del_node_attr(info_group, "numberOfPartitions")
            h5f.del_node_attr(info_group, "numberOfSubpartitions")

    if '/RGP' in h5f and (gene_families or partition or rgp):
        logging.getLogger().info("Erasing the formerly computer RGP...")
        pangenome.status["predictedRGP"] = "No"
        status_group._v_attrs.predictedRGP = False
        h5f.remove_node("/", "RGP")

        h5f.del_node_attr(info_group, "numberOfRGP")

    if '/spots' in h5f and (gene_families or partition or rgp or spots):
        logging.getLogger().info("Erasing the formerly computed spots...")
        pangenome.status["spots"] = "No"
        status_group._v_attrs.spots = False
        h5f.remove_node("/", "spots")

        h5f.del_node_attr(info_group, "numberOfSpots")

    if '/modules' in h5f and (gene_families or partition or modules):
        logging.getLogger().info("Erasing the formerly computed modules...")
        pangenome.status["modules"] = "No"
        status_group._v_attrs.modules = False
        h5f.remove_node("/", "modules")

        h5f.del_node_attr(info_group, "numberOfModules")
        h5f.del_node_attr(info_group, "numberOfFamiliesInModules")
        for info in ['CloudSpecInModules', 'PersistentSpecInModules', 'ShellSpecInModules', 'numberOfFamiliesInModules',
                     'StatOfFamiliesInModules']:
            if info in info_group._v_attrs._f_list():
                h5f.del_node_attr(info_group, info)

    h5f.close()


def write_pangenome(pangenome: Pangenome, filename, force: bool = False, disable_bar: bool = False):
    """
    Writes or updates a pangenome file

    :param pangenome: pangenome object
    :param filename: HDF5 file to save pangenome
    :param force: force to write on pangenome if information already exist
    :param disable_bar: Allow to disable progress bar
    """

    if pangenome.status["genomesAnnotated"] == "Computed":
        compression_filter = tables.Filters(complevel=1, shuffle=True, bitshuffle=True, complib='blosc:zstd')
        h5f = tables.open_file(filename, "w", filters=compression_filter)
        logging.getLogger().info("Writing genome annotations...")

        write_annotations(pangenome, h5f, disable_bar=disable_bar)

        pangenome.status["genomesAnnotated"] = "Loaded"
        h5f.close()
    elif pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"]:
        pass
    else:
        # if the pangenome is not Computed or not Loaded, it's probably not really in a good state
        # (or something new was coded).
        raise NotImplementedError("Something REALLY unexpected and unplanned for happened here. "
                                  "Please post an issue on github with what you did to reach this error.")

    # from there, appending to existing file.
    h5f = tables.open_file(filename, "a")

    if pangenome.status["geneSequences"] == "Computed":
        logging.getLogger().info("writing the protein coding gene dna sequences")
        write_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["geneSequences"] = "Loaded"

    if pangenome.status["genesClustered"] == "Computed":
        logging.getLogger().info("Writing gene families and gene associations...")
        write_gene_families(pangenome, h5f, force, disable_bar=disable_bar)
        logging.getLogger().info("Writing gene families information...")
        write_gene_fam_info(pangenome, h5f, force, disable_bar=disable_bar)
        if pangenome.status["genomesAnnotated"] in ["Loaded", "inFile"] and \
                pangenome.status["defragmented"] == "Computed":
            # if the annotations have not been computed in this run,
            # and there has been a clustering with defragmentation, then the annotations can be updated
            update_gene_fragments(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["genesClustered"] = "Loaded"
    if pangenome.status["neighborsGraph"] == "Computed":
        logging.getLogger().info("Writing the edges...")
        write_graph(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["neighborsGraph"] = "Loaded"
    if pangenome.status["partitioned"] == "Computed" and \
            pangenome.status["genesClustered"] in ["Loaded", "inFile"]:  # otherwise, it's been written already.
        update_gene_fam_partition(pangenome, h5f, disable_bar=disable_bar)
        pangenome.status["partitioned"] = "Loaded"

    if pangenome.status['predictedRGP'] == "Computed":
        logging.getLogger().info("Writing Regions of Genomic Plasticity...")
        write_rgp(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status['predictedRGP'] = "Loaded"

    if pangenome.status["spots"] == "Computed":
        logging.getLogger().info("Writing Spots of Insertion...")
        write_spots(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status['spots'] = "Loaded"

    if pangenome.status["modules"] == "Computed":
        logging.getLogger().info("Writing Modules...")
        write_modules(pangenome, h5f, force, disable_bar=disable_bar)
        pangenome.status["modules"] = "Loaded"

    write_status(pangenome, h5f)
    write_info(pangenome, h5f)

    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")
