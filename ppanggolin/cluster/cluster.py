#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
import logging
import tempfile
import subprocess
from collections import defaultdict
import os


#installed libraries
import networkx

#local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Gene
from ppanggolin.utils import read_compressed_or_not, getCurrentRAM
from ppanggolin.formats import writePangenome, readPangenome, writeCDSSequences, getStatus
from ppanggolin.annotate import genetic_codes, translate

def alignRep(faaFile, tmpdir, cpu):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir.name)#create a tmpdir in the tmpdir provided.
    seqdb =  tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","createdb",faaFile.name, seqdb.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    alndb =  tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","search",seqdb.name , seqdb.name, alndb.name, newtmpdir.name, "-a","--min-seq-id", "0.8", "-c", "0.8", "--cov-mode", "1", "--threads", str(cpu)]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Aligning cluster representatives...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    outfile =  tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","convertalis", seqdb.name ,seqdb.name, alndb.name, outfile.name,"--format-output","query,target,qlen,tlen,bits"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Extracting alignments...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    seqdb.close()
    alndb.close()
    newtmpdir.cleanup()
    return outfile

def firstClustering(sequences, tmpdir, cpu, code ):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir.name)#create a tmpdir in the tmpdir provided.
    seqNucdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","createdb"] 
    cmd.append(sequences)
    cmd.extend([seqNucdb.name,"--dont-shuffle","false"])
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Creating sequence database...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","translatenucs", seqNucdb.name, seqdb.name, "--threads", str(cpu), "--translation-table",code]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Translating sequences...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    cludb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","cluster",seqdb.name, cludb.name, newtmpdir.name , "--min-seq-id", "0.8", "-c", "0.8", "--threads", str(cpu), "--kmer-per-seq","80","--max-seqs","300"]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Clustering sequences...")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    logging.getLogger().info("Extracting cluster representatives...")
    repdb = tempfile.NamedTemporaryFile(mode="w", dir = newtmpdir.name)
    cmd = ["mmseqs","result2repseq", seqdb.name, cludb.name, repdb.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    reprfa = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","result2flat",seqdb.name, seqdb.name, repdb.name, reprfa.name]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    outtsv = tempfile.NamedTemporaryFile(mode="w", dir = tmpdir.name)
    cmd = ["mmseqs","createtsv",seqdb.name, seqdb.name, cludb.name,outtsv.name]
    logging.getLogger().debug(" ".join(cmd))
    logging.getLogger().info("Writing gene to family informations")
    subprocess.run(cmd, stdout=open(os.devnull,"w"))
    repdb.close()
    seqdb.close()
    cludb.close()
    seqNucdb.close()
    newtmpdir.cleanup()#deleting temporary directory.
    return reprfa, outtsv

def read_faa(faFileName):
    fam2seq = {}
    head = ""
    with open(faFileName.name,"r") as faFile:
        for line in faFile:
            if line.startswith('>'):
                head = line[1:].strip()
            else:
                fam2seq[head] = line.strip()
    return fam2seq

def read_tsv(tsvfileName):
    # reading tsv file
    genes2fam = {}
    fam2genes = defaultdict(set)
    with open(tsvfileName.name, "r") as tsvfile:
        for line in tsvfile:
            line = line.split()
            genes2fam[line[1]] = (line[0],None)#fam id, and its a gene (and not a fragment)
            fam2genes[line[0]].add(line[1])
    return genes2fam, fam2genes

def refineClustering(tsv, alnFile, fam2seq):
    simgraph = networkx.Graph()
    genes2fam, fam2genes = read_tsv(tsv)
    logging.getLogger().info(f"Starting with {len(fam2seq)} families")
    #create the nodes
    for fam, genes in fam2genes.items():
        simgraph.add_node(fam, nbgenes = len(genes) )
    #add the edges
    with open(alnFile.name,"r") as alnfile:
        for line in alnfile:
            line = line.split()
            
            if line[0] != line[1]:
                simgraph.add_edge(line[0],line[1], score = float(line[4]))
                simgraph.nodes[line[0]]["length"] = int(line[2])
                simgraph.nodes[line[1]]["length"] = int(line[3])
    for node, nodedata in simgraph.nodes(data = True):
        choice = (None, 0, 0, 0)
        for neighbor in simgraph.neighbors(node):
            nei = simgraph.nodes[neighbor]
            score = simgraph[neighbor][node]["score"]
            if nei["length"] > nodedata["length"] and nei["nbgenes"] >= nodedata["nbgenes"] and  choice[3] < score:
                choice = (genes2fam[neighbor][0], nei["length"] , nei["nbgenes"], score)#genes2fam[neighbor] instead of just neighbor in case that family has been assigned already (this is for smaller fragments that are closer to other fragments than the actual gene family)
        if choice[0] is not None:
            genestochange = fam2genes[node]
            for gene in genestochange:
                genes2fam[gene] = (choice[0], "F")
                fam2genes[choice[0]].add(gene)
            del fam2genes[node]
    newFam2seq = {}
    for fam in fam2genes:
        newFam2seq[fam] = fam2seq[fam]
    logging.getLogger().info(f"Ending with {len(newFam2seq)} gene families")
    return genes2fam, newFam2seq

def read_gene2fam(pangenome, gene2fam):
    logging.getLogger().info("Adding genes to the gene families")
    for gene, (family, is_frag) in gene2fam.items():
        fam = pangenome.addGeneFamily(family)
        geneObj = Gene(gene)
        geneObj.is_fragment = is_frag
        fam.addGene(geneObj)

def read_fam2seq(pangenome, fam2seq):
    logging.getLogger().info("Adding protein sequences to the gene families")
    for family, protein in fam2seq.items():
        fam = pangenome.addGeneFamily(family)
        fam.addSequence(protein)

def clustering(pangenome, sequences, tmpdir, cpu , nodefrag = False, code = "11"):
    newtmpdir = tempfile.TemporaryDirectory(dir = tmpdir)
    logging.getLogger().info("Clustering all of the genes sequences...")
    rep, tsv = firstClustering(sequences, newtmpdir, cpu, code)
    fam2seq = read_faa(rep)
    if nodefrag:
        genes2fam = read_tsv(tsv)[0]
        tsv.close()
        rep.close()
        newtmpdir.cleanup()
    else:
        logging.getLogger().info("Associating fragments to their original gene family...")
        aln = alignRep(rep, newtmpdir, cpu)
        genes2fam, fam2seq = refineClustering(tsv, aln, fam2seq)
        aln.close()
        tsv.close()
        rep.close()
        newtmpdir.cleanup()
    read_fam2seq(pangenome, fam2seq)
    read_gene2fam(pangenome, genes2fam)
    pangenome.status["genesClustered"] = "Computed"
    pangenome.status["geneFamilySequences"] = "Computed"

def launch(args):
    """ launch the clustering step"""

    # pangenome = readPangenome(args.pangenome, cpu = args.cpu, annotation = True)#load the annotations
    # pangenome.stats()
    logging.getLogger().debug(f"Ram used at the start : {getCurrentRAM()}")
    tmpFile = tempfile.NamedTemporaryFile(mode="w", dir = args.tmpdir)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    writeCDSSequences(pangenome, tmpFile, args.cpu)
    logging.getLogger().info("Done extracting all the protein sequences")
    clustering(pangenome, tmpFile.name, args.tmpdir, args.cpu, args.nodefrag, args.translation_table)
    writePangenome(pangenome, pangenome.file)
    tmpFile.close()
    logging.getLogger().debug(f"RAM used after clustering : {getCurrentRAM()}")
    ##now deal with those things (like write them ?)
    logging.getLogger().info("Done with the clustering")

def clusterSubparser(subparser):
    parser = subparser.add_parser("cluster",help = "Cluster proteins in protein families")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--nodefrag', required=False,default=False, action="store_true", help = "Do not use the defragmentation strategy to associated potential fragments with their original gene family.")
    optional.add_argument('--tsv', required=False, type=str, help="A tab-separated file listing the gene cluster IDs (which must be the representative gene ID), and the gene IDs in the cluster. Optionally, a 'F' can be added in a third column to indicate that the gene is a fragment. One gene per line.")
    optional.add_argument("--translation_table",required=False, default="11", help = "Translation table (genetic code) to use.")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('-p','--pangenome',  required=True, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism.")
    return parser