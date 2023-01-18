#!/usr/bin/env python3
# coding: utf8

# default libraries
import logging
from collections.abc import Iterable

# local libraries
from ppanggolin.genome import Organism
from ppanggolin.region import Region
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.edge import Edge
from ppanggolin.utils import canonical, get_file_size

class Node:
    """A class representing a VG node"""
    def __init__(self, sequence: str):
        self.len_sequence = len(sequence)   # we need to remind the length of the sequence of each node for statistical computations
        self.traversed_path = []            # ids (int) of the traversed paths. This corresponds to indexes in the Pangenome.paths list
        self.sequence = sequence
class Path:
    """A class representing a VG path composed of several nodes"""
    def __init__(self, ID):
        self.ID = ID
        self.node_ids   = []                # a set of node ids (ints)
        self.family = None              # id (string) of the cluster this path belongs to
        self.genes = set()                # ids of the genes that generated this path with their counts.
        self.organisms = {}                # dict of the organism that shared this path with their counts if it is present several times

    def get_sequence_length(self):
        return sum([self.nodes[node_id].len_sequence for node_id in self.node_ids])

class Pangenome:
    """
    This is a class representing your pan. It is used as a basic unit for all the analysis to access to the
    different elements of your pan, such as organisms, contigs, genes or gene families. It has setter and getter
    methods for most elements in your pan, and you can use those to add new elements to it,
    or get objects that have a specific identifier to manipulate them directly.
    """

    def __init__(self):
        """Constructor method.
        """
        self.file = None

        # basic parameters
        self._famGetter = {}
        self._org_index = None
        self._fam_index = None
        self.max_fam_id = 0
        self._orgGetter = {}
        self._edgeGetter = {}
        self._regionGetter = {}
        self.spots = set()
        self.modules = set()
        ###
        self.nodes = {}                     # key: id (unsigned int), value: Node 
        self.paths = []                     # a list of ordered paths (id of a path is its ranks in this list)
        self.paths_name_to_ids = {}         # no choice: each path has a string name, eg gi|1388876906|ref|NZ_CP028116.1|_1000. Two identical paths (eg 684619+,684620+,684618+) may have distinct occurrences and thus names (as comming from distinct genes). Hence one storesfor each path name its unique path id.
        self.paths_content_to_ids = {}      # no choice: each path has a content, eg 684619+,684620+,684618+. This is the key to know its UNIQUE id (int) that is also the rank in the self.paths list
        self.species_names = set()          # store all species ids NZ_CP007592.1, NC_013654.1, ...

        self.status = {
            'genomesAnnotated': "No",
            'geneSequences': "No",
            'genesClustered': "No",
            'defragmented': "No",
            'geneFamilySequences': "No",
            'neighborsGraph': "No",
            'partitioned': "No",
            'predictedRGP': "No",
            'spots': "No",
            'modules': 'No',
            'variation_graphs': 'No'
        }
        self.parameters = {}

    def add_file(self, pangenome_file):
        """Links an HDF5 file to the pan. If needed elements will be loaded from this file,
        and anything that is computed will be saved to this file when
        :func:`ppanggolin.formats.writeBinaries.writePangenome` is called.

        :param pangenome_file: A string representing the filepath to the hdf5 pan file
        to be either used or created
        :type pangenome_file: str
        """
        from ppanggolin.formats.readBinaries import get_status
        # importing on call instead of importing on top to avoid cross-reference problems.
        get_status(self, pangenome_file)
        self.file = pangenome_file

    """ Gene Methods"""
    @property
    def genes(self):
        """Creates the geneGetter if it does not exist, and returns all the genes of all organisms in the pangenome.
        
        :return: list of :class:`ppanggolin.genome.Gene`
        :rtype: list
        """
        try:
            return list(self._geneGetter.values())
        except AttributeError:  # in that case the gene getter has not been computed
            self._mk_gene_getter()  # make it
            return self.genes  # return what was expected

    def _yield_genes(self):
        """
        Use a generator to get all the genes of a pan

        :return: an iterator of :class:`ppanggolin.genome.Gene`
        :rtype: Iterator[:class:`ppanggolin.genome.Gene`]
        """
        if self.number_of_organisms() > 0:  # if we have organisms, they're supposed to have genes
            for org in self.organisms:
                for contig in org.contigs:
                    for gene in contig.genes:
                        yield gene
        elif self.number_of_gene_families() > 0:
            # we might have no organism loaded, in that case there are gene families.
            for geneFam in self.gene_families:
                for gene in geneFam.genes:
                    yield gene

    def _mk_gene_getter(self):
        """
            Builds the :attr:`ppanggolin.pan.Pangenome._geneGetter` of the pan

            Since the genes are never explicitly 'added' to a pan (but rather to a gene family, or a contig),
            the pan cannot directly extract a gene from a geneID since it does not 'know' them.
            if at some point we want to extract genes from a pan we'll create a geneGetter.
            The assumption behind this is that the pan has been filled and no more gene will be added.
        """
        self._geneGetter = {}
        for gene in self._yield_genes():
            self._geneGetter[gene.ID] = gene

    def get_gene(self, gene_id):
        """returns the gene that has the given `geneID`

        :param gene_id: The gene ID to look for
        :type gene_id: any
        :return: returns the gene that has the ID `geneID`
        :rtype: :class:`ppanggolin.genome.Gene`
        :raises KeyError: If the `geneID` is not in the pan
        """
        try:
            return self._geneGetter[gene_id]
        except AttributeError:
            # in that case, either the gene getter has not been computed, or the geneID is not in the pan.
            self._mk_gene_getter()  # make it
            return self.get_gene(
                gene_id)  # return what was expected. If the geneID does not exist it will raise an error.
        except KeyError:
            raise KeyError(f"{gene_id} does not exist in the pan.")

    """Gene families methods"""
    @property
    def gene_families(self):
        """returns all the gene families in the pan
        
        :return: list of :class:`ppanggolin.geneFamily.GeneFamily`
        :rtype: list
        """
        return list(self._famGetter.values())

    def _create_gene_family(self, name):
        """Creates a gene family object with the given `name`

        :param name: the name to give to the gene family. Must not exist already.
        :type name: any
        :return: the created GeneFamily object
        :rtype: :class:`ppanggolin.geneFamily.GeneFamily`
        """
        new_fam = GeneFamily(family_id=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._famGetter[new_fam.name] = new_fam
        return new_fam

    def number_of_gene_families(self):
        """Returns the number of gene families present in the pan

        :return: the number of gene families
        :rtype: int
        """
        return len(self._famGetter)

    def get_gene_family(self, name):
        """returns the gene family that has the given `name`

        :param name: The gene family name to look for
        :type name: any
        :return: returns the gene family that has the name `name`
        :rtype: :class:`ppanggolin.geneFamily.GeneFamily`
        """
        return self._famGetter[name]

    def add_gene_family(self, name):
        """
            Get the :class:`ppanggolin.geneFamily.GeneFamily` object that has the given `name`. If it does not exist,
            creates it.
            returns the geneFamily object.

            :param name: The gene family name to get if it exists, and create otherwise.
            :type name: str
        """
        fam = self._famGetter.get(name)
        if fam is None:
            fam = self._create_gene_family(name)
        return fam

    """Graph methods"""
    @property
    def edges(self):
        """returns all the edges in the pan graph
        
        :return: list of :class:`ppanggolin.pan.Edge`
        :rtype: list
        """
        return list(self._edgeGetter.values())

    def add_edge(self, gene1, gene2):
        """
        Adds an edge between the two gene families that the two given genes belong to. Genes object are expected,
        and they are also expected to have a family assigned

        :param gene1: The first gene
        :type gene1: :class:`ppanggolin.genome.Gene`
        :param gene2: The second gene
        :type gene2: :class:`ppanggolin.genome.Gene`
        :return: the created Edge
        :rtype: :class:`ppanggolin.pangenome.Edge`
        """
        key = frozenset([gene1.family, gene2.family])
        edge = self._edgeGetter.get(key)
        if edge is None:
            edge = Edge(gene1, gene2)
            self._edgeGetter[key] = edge
        else:
            edge.add_genes(gene1, gene2)
        return edge

    """Organism methods"""
    @property
    def organisms(self):
        """returns all the organisms in the pan
        
        :return: list of :class:`ppanggolin.genome.Organism`
        :rtype: list
        """
        return list(self._orgGetter.values())

    def number_of_organisms(self):
        """Returns the number of organisms present in the pan
        
        :return: the number of organism
        :rtype: int
        """
        return len(self._orgGetter)

    def get_organism(self, org_name):
        """
        Get an organism that is expected to be in the pan using its name, which is supposedly unique.
        Raises an error if the organism does not exist.

        :param org_name: Name of the :class:`ppanggolin.genome.Organism` to get
        :type org_name: str
        :return: The related Organism object
        :rtype: :class:`ppanggolin.genome.Organism`
        :raises KeyError: If the provided name is not in the pan
        """
        try:
            return self._orgGetter[org_name]
        except KeyError:
            raise KeyError(f"{org_name} does not seem to be in your pan")

    def add_organism(self, new_org):
        """
        adds an organism that did not exist previously in the pan if an :class:`ppanggolin.genome.Organism`
        object is provided. If an organism with the same name exists it will raise an error.
        If a :class:`str` object is provided, will return the corresponding organism that has this name
        OR create a new one if it does not exist.

        :param new_org: Organism to add to the pan
        :type new_org: :class:`ppanggolin.genome.Organism` or str
        :return: The created organism
        :rtype: :class:`ppanggolin.genome.Organism`
        :raises TypeError: if the provided `newOrg` is neither a str nor a :class:`ppanggolin.genome.Organism`
        """
        if isinstance(new_org, Organism):
            old_len = len(self._orgGetter)
            self._orgGetter[new_org.name] = new_org
            if len(self._orgGetter) == old_len:
                raise KeyError(f"Redondant organism name was found ({new_org.name})."
                               f"All of your organisms must have unique names.")
        elif isinstance(new_org, str):
            org = self._orgGetter.get(new_org)
            if org is None:
                org = Organism(new_org)
                self._orgGetter[org.name] = org
            new_org = org
        else:
            raise TypeError("Provide an Organism object or a str that will serve as organism name")
        return new_org

    def get_index(self):  # will not make a new index if it exists already
        """Creates an index for Organisms (each organism is assigned an Integer).

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if self._org_index is None:  # then the index does not exist yet
            self._org_index = {}
            for index, org in enumerate(self.organisms):
                self._org_index[org] = index
        return self._org_index

    def compute_family_bitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pan.Pangenome.getIndex`, generated a bitarray
        for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionnary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if self._org_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_index()
        for fam in self.gene_families:
            fam.mk_bitarray(self._org_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._org_index

    def get_fam_index(self):  # will not make a new index if it exists already
        """Creates an index for gene families (each family is assigned an Integer).

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if self._fam_index is None:  # then the index does not exist yet
            self._fam_index = {}
            for index, fam in enumerate(self.gene_families):
                self._fam_index[fam] = index
        return self._fam_index

    def compute_org_bitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pan.Pangenome.get_fam_index`, generated a bitarray
        for each gene family.
        If the family j is present in the organism with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if self._fam_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for org in self.organisms:
            org.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index

    """RGP methods"""
    @property
    def regions(self):
        """returns all the regions (RGP) in the pan

        :return: list of :class:`ppanggolin.region.Region`
        :rtype: list
        """
        return list(self._regionGetter.values())

    def get_or_add_region(self, region_name):
        """Returns a region with the given `regionName`. Creates it if it does not exist.

        :param region_name: The name of the region to return
        :type region_name: str
        :return: The region
        :rtype: :class:`ppanggolin.region.Region`
        """
        try:
            return self._regionGetter[region_name]
        except KeyError:  # then the region is not stored in this pan.
            new_region = Region(region_name)
            self._regionGetter[region_name] = new_region
            return new_region

    def get_multigenics(self, dup_margin, persistent=True):
        """
        Returns the multigenic persistent families of the pan graph. A family will be considered multigenic
        if it is duplicated in more than `dup_margin` of the genomes where it is present.

        :param dup_margin: the ratio of presence in multicopy above which a gene family is considered multigenic
        :type dup_margin: float
        :param persistent: if we consider only the persistent genes
        :type persistent: bool
        :return: a `set` of gene families considered multigenic
        :rtype: set[:class:`ppanggolin.geneFamily.GeneFamily`]
        """
        multigenics = set()
        for fam in self.gene_families:
            if fam.named_partition == "persistent" or not persistent:
                dup = len([genes for org, genes in fam.get_org_dict().items() if
                           len([gene for gene in genes if not gene.is_fragment]) > 1])
                if (dup / len(fam.organisms)) >= dup_margin:  # tot / nborgs >= 1.05
                    multigenics.add(fam)
        # logging.getLogger().info(f"{len(multigenics)} gene families are defined as being multigenic.
        # (duplicated in more than {dup_margin} of the genomes)")
        return multigenics

    def add_regions(self, region_group):
        """Takes an Iterable or a Region object and adds it to the pan

        :param region_group: a region or an Iterable of regions to add to the pan
        :type region_group: :class:`ppanggolin.region.Region` or Iterable[:class:`ppanggolin.region.Region`]
        :raises TypeError: if regionGroup is neither a Region nor an Iterable[:class:`ppanggolin.region.Region`]
        """
        old_len = len(self._regionGetter)
        if isinstance(region_group, Iterable):
            for region in region_group:
                self._regionGetter[region.name] = region
            if len(self._regionGetter) != len(region_group) + old_len:
                raise Exception("Two regions had an identical name, which was unexpected.")
        elif isinstance(region_group, Region):
            self._regionGetter[region_group.name] = region_group
        else:
            raise TypeError(f"An iterable or a 'Region' type object were expected, "
                            f"but you provided a {type(region_group)} type object")

    """Spot methods"""
    def add_spots(self, spots):
        """Adds the given iterable of spots to the pan.

        :param spots: An iterable of :class:`ppanggolin.region.Spot`.
        :type spots: Iterable[:class:`ppanggolin.region.Spot`]
        """
        self.spots |= set(spots)

    """Modules methods"""
    def add_modules(self, modules):
        """Adds the given iterable of modules to the pan

        :param modules: an iterable of :class:`ppanggolin.module.Module`
        :type modules: Iterable[:class:`ppanggolin.module.Module`]
        """
        self.modules |= set(modules)

    def compute_mod_bitarrays(self, part='all'):
        """Based on the index generated by :meth:`ppanggolin.pan.Pangenome.get_fam_index`, generated a bitarray
        for each gene family.
        If the family j is present in the module with the index i, the bit at position i will be 1. If it is not,
        the bit will be 0.
        The bitarrays are gmpy2.xmpz object.

        :return: A dictionary with :class:`ppanggolin.genome.Organism` as key and `int` as value.
        :rtype: dict[:class:`ppanggolin.genome.Organism`, int]
        """
        if self._fam_index is None:
            # then the bitarrays don't exist yet, since the org index does not exist either.
            self.get_fam_index()
        for module in self.modules:
            module.mk_bitarray(index=self._fam_index, partition=part)
        # case where there is an index but the bitarrays have not been computed???
        return self._fam_index

    def fill_from_GFA_file(self, gfa_file_name: str):
        """
        PARSE GFA FILE
        From the gfa file, create a dictionary of the VG nodes and paths in the graph.
        S lines contain node ID and its sequence
        P lines contain path ID, list of node ID with orientation and cover of the nodes
        L lines contain links between nodes (we dont care)
        """
        logging.getLogger().debug("Load the GFA file")
        path_id = 0
        with open(gfa_file_name, 'r') as gfa_file:
            size_file = get_file_size(gfa_file)
            while True:
                line = gfa_file.readline()

                # end of the file
                if not line:
                    break

                current_seek = gfa_file.tell()
      
                line = line.strip().split('\t')
                # if line S, create node
                if line[0] == 'S':
                    # S       1       ACCACGATTACGCTGGCGCTTA
                    self.nodes[int(line[1])] = Node(line[2]) 
                # if line P, create paths and add paths infos in nodes
                elif line[0] == 'P':
                    # P       gi|1388876906|ref|NZ_CP028116.1|_1000   684619+,684620+,684618+ 187M,187M,1M 
                    # path_id = line[1]
                    str_node_list = canonical(line[2])
                    gene = self.get_gene(line[1])
                    # If this path was already seen, we simply add this strain_id to the path.strain_ids
                    if str_node_list in self.paths_content_to_ids:
                        already_seen_path_id = self.paths_content_to_ids[str_node_list]
                        if gene not in self.paths[already_seen_path_id].genes:
                            self.paths[already_seen_path_id].organisms[gene.organism]=0
                        self.paths[already_seen_path_id].organisms[gene.organism]+=1
                        self.paths[already_seen_path_id].genes.add(gene)
                        self.paths_name_to_ids[line[1]] = already_seen_path_id
                        continue # nothing more to do, no incrementation of path_id, as no new path was created

                    path = Path(path_id)
                    path.family=gene.family
                    if gene.organism not in path.organisms:
                        path.organisms[gene.organism]=0
                    path.organisms[gene.organism]+=1
                    node_list = str_node_list.split(',')
                    # seq = ""
                    for node_info in node_list:
                        node_id = int(node_info[:-1])
                        node = self.nodes[node_id]
                        # add the current path id to the traversed paths of this node
                        # we could have used a set for traversed_path, but this requires more memory
                        if path_id not in node.traversed_path:
                            node.traversed_path.append(path_id)

                        path.node_ids.append(node_id)
                        # path.hamming_distance.append(0)
                    # assert line[1] not in self.paths_name_to_ids # todo: remove when tested that this path ids has never been seen before

                    self.paths.append(path)                                 # store this new path
                    self.paths_name_to_ids[line[1]] = path_id
                    self.paths_content_to_ids[str_node_list] = path_id
                    path_id+=1
        self.status["variation_graphs"]="Loaded"

    def get_matching_path(self, path_as_node_list):
        """
        INPUT = a node list (from alignement)
        to speed the search we get the starting node of the alignment and search only on paths crossing this start node
        to speed the seach, we use the node position in the selected path and extend the list of nodes by the same length of the alignment
        there is a match if the seed+extend list is identical to the node list of the alignment
        OUPUT = list of paths where the node list match
        """
        if self.status["variation_graphs"] != "No":
            result = []  # (path_id, corresponding starting node_id, number of nodes mapped)
            start_node = path_as_node_list[0]
            paths_sel = self.nodes[start_node].traversed_path

            # for each path crossing the start node
            for path_id in paths_sel:
                colored_path = self.paths[path_id]
                all_pos_first = [i for i, v in enumerate(colored_path.node_ids) if
                                 v == start_node]  # start_node from alignment may have several occurrences in the colored path

                for pos_first in all_pos_first:
                    if colored_path.node_ids[pos_first:pos_first + len(path_as_node_list)] == path_as_node_list:
                        result.append((path_id, pos_first))

            return result

    def get_sequence_length(self, path):
        return sum([self.nodes[node_id].len_sequence for node_id in path.node_ids])