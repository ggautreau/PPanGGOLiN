#! /usr/bin/env python3

import pytest
from random import choices, randint, sample

from ppanggolin.genome import Gene, Organism
from ppanggolin.pangenome import Edge, Pangenome
from ppanggolin.geneFamily import GeneFamily


def test_cstr():
    o_pang = Pangenome()
    assert isinstance(o_pang, Pangenome)

    for attr in "max_fam_id", "parameters", "status":
        assert hasattr(o_pang, attr)
    assert o_pang.max_fam_id == 0
    assert o_pang.parameters == {}
    assert o_pang.status == {
        'genomesAnnotated': "No",
        'geneSequences': "No",
        'genesClustered': "No",
        'defragmented': "No",
        'geneFamilySequences': "No",
        'neighborsGraph': "No",
        'partitioned': "No",
        'predictedRGP': "No",
        'spots': "No",
        'modules': "No"
    }


@pytest.fixture
def o_pang():
    return Pangenome()


# @pytest.mark.xfail(reason="not implemented !")
# def test_add_file(o_pang):
#     assert False  # need to generate a valid file several time

@pytest.fixture
def l_orgs():
    l_orgs = []
    for i_org in range(randint(5, 20)):
        o_org = Organism(str(i_org))
        l_orgs.append(o_org)

    return l_orgs


def test_organisms(o_pang, l_orgs):
    # 'set' because order is not guaranted
    # and org should be unique
    assert set(o_pang.organisms) == set()

    # add Org from Org
    for o_org in l_orgs:
        o_pang.add_organism(o_org)

    # add Org from string
    for i_org in range(randint(5, 20)):
        o_org = o_pang.add_organism(str(i_org))
        l_orgs.append(o_org)

    assert set(o_pang.organisms) == set(l_orgs)


def test_add_organism_str(o_pang):
    o_org = o_pang.add_organism("org1")
    assert o_org in o_pang.organisms
    assert isinstance(o_org, Organism)
    assert set(o_pang.organisms) == {o_org}


def test_add_organism(o_pang):
    o_org = Organism("org")
    assert o_pang.add_organism(o_org) == o_org
    assert set(o_pang.organisms) == {o_org}


def test_number_of_organism(o_pang, l_orgs):
    assert o_pang.number_of_organisms() == 0

    for o_org in l_orgs:
        o_pang.add_organism(o_org)

    assert o_pang.number_of_organisms() == len(l_orgs)


def test_add_gene_family_one(o_pang):
    name = "fam1"
    o_fam1 = o_pang.add_gene_family(name)
    assert isinstance(o_fam1, GeneFamily)
    assert 1 == o_pang.max_fam_id


def test_add_gene_family_same(o_pang):
    name = "fam1"
    o_fam1 = o_pang.add_gene_family(name)
    o_fam2 = o_pang.add_gene_family(name)
    assert o_fam1 == o_fam2


def test_add_gene_family_many(o_pang):
    n_fams = randint(5, 20)
    for i_fam in range(n_fams):
        o_pang.add_gene_family(str(i_fam))
    assert n_fams == o_pang.max_fam_id


def test_get_gene_family(o_pang):
    name = "fam1"
    o_fam = o_pang.add_gene_family(name)
    assert o_pang.get_gene_family(name) == o_fam

    for i_fam in range(randint(5, 20)):
        o_pang.add_gene_family(str(i_fam))
    # still true after many insert
    assert o_pang.get_gene_family(name) == o_fam


def test_number_of_gene_families_empty(o_pang):
    assert o_pang.number_of_gene_families() == 0


def test_number_of_gene_families(o_pang):
    n_fams = randint(5, 10)

    for i_fam in sample(range(20), k=n_fams):
        o_pang.add_gene_family(str(i_fam))
    assert o_pang.number_of_gene_families() == n_fams


def test_gene_families_empty(o_pang):
    # 'set' because order is not guaranted
    assert set(o_pang.gene_families) == set()


def test_gene_families(o_pang):
    l_ints = choices(range(20), k=10)
    s_fams = set()
    for i_fam in l_ints:
        o_fam = o_pang.add_gene_family(str(i_fam))
        s_fams.add(o_fam)

    assert set(o_pang.gene_families) == s_fams


def test_genes_empty(o_pang):
    assert list(o_pang.genes) == []


# code copy-pasted from test_Edge.py
@pytest.fixture()
def make_gene_pair():
    def _make_gene_pair(org, gene_id1, gene_id2):
        """create a pair of genes that belong to the same organism."""
        lo_genes = []
        for k in gene_id1, gene_id2:
            o_gene = Gene(k)
            o_gene.fill_parents(org, None)

            lo_genes.append(o_gene)

            o_family = GeneFamily(k, k)
            o_family.add_gene(o_gene)
        return tuple(lo_genes)
    return _make_gene_pair


@pytest.fixture()
def make_org_with_genes():
    def _make_org_with_genes(org):
        """make an organism, add from 2 to 10 contigs
        with 2 to 10 genes each."""
        l_genes = []
        o_org = Organism(org)
        for i in range(randint(2, 10)):
            o_ctg = o_org.get_contig("k_{}".format(i))
            for j in range(randint(2, 10)):
                name = "{}.{}.{}".format(org, o_ctg.name, j)
                o_gene = Gene(name)
                o_gene.position = j
                o_gene.start = j
                o_ctg.add_gene(o_gene)
                l_genes.append(o_gene)
        return o_org, l_genes

    return _make_org_with_genes


@pytest.fixture()
def fill_fam_with_genes():
    def _fill_fam_with_genes(o_fam):
        """add genes with names from 2 to 10 to a geneFamily object."""
        l_genes = []
        for i in range(2, 10):
            name = "{}_{}".format(o_fam.name, i)
            o_gene = Gene(name)
            o_fam.add_gene(o_gene)
            l_genes.append(o_gene)
        return l_genes

    return _fill_fam_with_genes


def test_genes_organism_debug(o_pang, make_org_with_genes):
    # orgs with genes.
    o_org, l_genes = make_org_with_genes("org1")
    o_pang.add_organism(o_org)
    l_expected = sorted(l_genes, key=lambda g: g.ID)
    l_observed = sorted(o_pang.genes, key=lambda g: g.ID)
    assert l_observed == l_expected


def test_genes_genefamilies(o_pang, fill_fam_with_genes):
    """Genes are added in pan through their family."""
    # geneFamily with genes.
    o_fam = o_pang.add_gene_family("fam1")
    l_genes = fill_fam_with_genes(o_fam)  # the list of genes, and the geneFam are supposed to be the same
    l_expected = sorted(l_genes, key=lambda g: g.ID)
    l_observed = sorted(o_pang.genes, key=lambda g: g.ID)
    print(o_pang.genes)
    assert l_observed == l_expected


def test_edges_empty(o_pang):
    assert list(o_pang.edges) == []


def test_add_edge(o_pang, make_gene_pair):
    name = "gene_fam"  # gene/fam name
    to_genes = make_gene_pair("org", name, name)

    o_edge1 = o_pang.add_edge(*to_genes)
    assert isinstance(o_edge1, Edge)

    # addEdge doesn't act the same when the edge already exists.
    o_edge2 = o_pang.add_edge(*to_genes)
    assert o_edge2 == o_edge1


def test_edges_one(o_pang, make_gene_pair):
    name = "gene_fam"  # gene/fam name
    to_genes = make_gene_pair("org", name, name)

    lo_edges = []
    n = randint(1, 5)
    for _ in range(n):
        lo_edges.append(o_pang.add_edge(*to_genes))

    # always the same family couple
    #    = one edge, with several couple of genes
    # I use set because edges are uniques, it is not a multigraph.
    assert set(o_pang.edges) == set(lo_edges)
    assert len(o_pang.edges) == 1

    o_edge = list(o_pang.edges).pop()
    assert o_edge.gene_pairs == [to_genes for _ in range(n)]


def test_edges_many_rand(o_pang, make_gene_pair):
    lo_edges = []
    n = randint(1, 5)
    for i in range(n):
        name1 = "gene_" + str(i)  # gene/fam name
        name2 = str(i) + "_gene"  # gene/fam name
        to_genes = make_gene_pair("org", name1, name2)
        lo_edges.append(o_pang.add_edge(*to_genes))
    # I use set because edges are uniques, it is not a supergraph.
    assert set(o_pang.edges) == set(lo_edges)


def test_edges_several(o_pang, make_gene_pair):
    # little more sophisticated
    to_genes = make_gene_pair("org", "g1", "g2")
    o_fam2 = to_genes[1].family
    o_pang.add_edge(*to_genes)

    to_genes = make_gene_pair("org", "g1", "g3")
    o_fam3 = to_genes[1].family
    o_pang.add_edge(*to_genes)
    # g3 -- g1 -- g2

    to_genes = make_gene_pair("org", "g22", "g33")
    o_fam2.add_gene(to_genes[0])
    o_fam3.add_gene(to_genes[1])
    o_pang.add_edge(*to_genes)
    # g2 -- g3

    assert len(o_pang.edges) == 3


def test_get_index(o_pang, l_orgs):
    for o_org in l_orgs:
        o_pang.add_organism(o_org)
    idx = o_pang.get_org_index()

    # after the method, the index exist
    assert o_pang.get_org_index() is idx

    # all orgs are in the index
    l_observed = sorted(idx.keys(), key=lambda x: x.name)
    l_orgs.sort(key=lambda x: x.name)
    assert l_observed == l_orgs


def test_compute_family_bitarrays(o_pang, l_orgs):
    for o_org in l_orgs:
        o_pang.add_organism(o_org)
    idx = o_pang.get_org_index()
    assert o_pang.compute_family_bitarrays() is idx


def test_family_have_bitarrays(o_pang, l_orgs):
    """test that after the method all the families have a bitarray."""
    n_fams = randint(5, 10)

    l_fams = []
    for i_fam in sample(range(20), k=n_fams):
        l_fams.append(o_pang.add_gene_family(str(i_fam)))
    o_pang.compute_family_bitarrays()
    for o_fam in l_fams:
        assert hasattr(o_fam, 'bitarray')


def test_get_gene_empty(o_pang):
    with pytest.raises(KeyError):
        o_pang.get_gene(33)


def test_get_gene_org(o_pang, make_org_with_genes):
    # orgs with genes.
    o_org, l_genes = make_org_with_genes("org")
    o_pang.add_organism(o_org)

    n = len(l_genes)
    for o_gene in sample(l_genes, randint(4, n)):
        assert o_pang.get_gene(o_gene.ID) == o_gene


def test_get_gene_fam(o_pang, fill_fam_with_genes):
    o_fam = o_pang.add_gene_family("fam")
    l_genes = fill_fam_with_genes(o_fam)

    for o_gene in l_genes:
        assert o_pang.get_gene(o_gene.ID) == o_gene
