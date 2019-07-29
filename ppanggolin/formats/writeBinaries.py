#!/usr/bin/env python3
#coding:utf-8

#default libraries
from pathlib import Path
import logging
import os
import warnings

#installed libraries
from tqdm import tqdm
import tables

#local libraries
from ppanggolin.pangenome import Pangenome, GeneFamily, Edge
from ppanggolin.genome import Organism, Contig, Gene
from ppanggolin.utils import getCurrentRAM

def geneDesc(geneLength):
    return {
            'organism':tables.StringCol(itemsize=1024),
            "contig":{
                    'name':tables.StringCol(itemsize=1024),
                    "is_circular":tables.BoolCol()
            },
            "gene":{ 
                'ID':tables.StringCol(itemsize=128),
                'start':tables.UInt64Col(),
                'stop':tables.UInt64Col(),
                'strand':tables.StringCol(itemsize=1),
                'type':tables.StringCol(itemsize=32),
                'position':tables.UInt32Col(),
                'name':tables.StringCol(itemsize=1024),
                'product':tables.StringCol(itemsize=1024),
                'genetic_code':tables.UInt32Col(),
                'is_fragment':tables.BoolCol(dflt = False),
                'dna':tables.StringCol(itemsize=geneLength),#The longest know protein to date is something around 36 000 amino acids.
            }
            }

def writeAnnotations(pangenome, h5f):
    """
        Function writing all of the pangenome's annotations
    """
    # original_warnings = list(warnings.filters)
    # warnings.simplefilter('ignore', tables.NaturalNameWarning)
    annotation = h5f.create_group("/","annotations","Annotations of the pangenome's organisms")
    nbGenes = 0
    maxlen = 0
    for gene in pangenome.genes:
        if len(gene.dna)>maxlen:
            maxlen = len(gene.dna)
        nbGenes+=1
    table = h5f.create_table(annotation, "genes", geneDesc(maxlen), expectedrows=nbGenes)
    bar = tqdm(pangenome.organisms, unit="genome")
    geneRow = table.row
    for org in bar:
        for contig in org.contigs:
            for gene in contig.genes:
                geneRow["organism"] = org.name
                geneRow["contig/name"] = contig.name
                geneRow["contig/is_circular"] = contig.is_circular#this should be somewhere else.
                geneRow["gene/ID"]= gene.ID
                geneRow["gene/start"] = gene.start
                geneRow["gene/stop"] = gene.stop
                geneRow["gene/strand"] = gene.strand
                geneRow["gene/type"] = gene.type
                geneRow["gene/position"] = gene.position
                geneRow["gene/name"] = gene.name
                geneRow["gene/product"] = gene.product
                geneRow["gene/is_fragment"] = gene.is_fragment
                geneRow["gene/genetic_code"] = gene.genetic_code
                if hasattr(gene, "dna"):
                    geneRow["gene/dna"] = gene.dna
                geneRow.append()
    table.flush()

    bar.close()
    # warnings.filters = original_warnings

def writeGeneFamilies(pangenome, h5f):
    """
        Function writing all of the pangenome's gene families
    """
    geneFamilies = h5f.create_group("/", "geneFamilies","Gene families of the pangenome")
    for geneFam in pangenome.geneFamilies:
        geneFamGroup = h5f.create_group(geneFamilies,geneFam.name)
        genesTable = h5f.create_table(geneFamGroup, "genes", {"organism":tables.StringCol(), "gene":tables.StringCol()})
        geneRow = genesTable.row
        for gene in geneFam.genes:
            geneRow["gene"] = gene.ID
            geneRow["organism"] = gene.organism.name
            geneRow.append()
        genesTable.flush()

def writeGraph(pangenome, h5f):
    edges = h5f.create_group("/", "graph","Edges of the neighbors graph")
    counter = 0#for naming uniquely edges, as they don't really have a proper way of identifying themselves.
    for edge in pangenome.edges:
        edgeGroup = h5f.create_group(edges, str(counter))
        edgeGroup.source = edge.source.name
        edgeGroup.target = edge.target.name
        edgeTable = h5f.create_table(edgeGroup, 'genes', { 'organism':tables.StringCol(), 'geneTarget':tables.StringCol(), 'geneSource':tables.StringCol() })
        counter +=1
        edgeRow = edgeTable.row
        for org, genePairs in edge.organisms.items():
            for gene1, gene2 in genePairs:
                edgeRow["organism"] = org.name
                edgeRow["geneTarget"] = gene1.name
                edgeRow["geneSource"] = gene2.name
                edgeRow.append()
        edgeTable.flush()

def writeStatus(pangenome, h5f):
    statusGroup = h5f.create_group("/","status","Statuses of the pangenome's content")
    statusGroup._v_attrs.genomesAnnotated = True if pangenome.status["genomesAnnotated"] in ["Yes","Loaded"] else False
    statusGroup._v_attrs.genesClustered = True if pangenome.status["genesClustered"] in ["Yes","Loaded"] else False
    statusGroup._v_attrs.NeighborsGraph = True if pangenome.status["NeighborsGraph"] in ["Yes","Loaded"] else False
    statusGroup._v_attrs.Partitionned = True if pangenome.status["Partitionned"] in ["Yes","Loaded"] else False

def writePangenome(pangenome, filename):
    """
        Writes or updates a pangenome file
        pangenome is the corresponding pangenome object, filename the h5 file and status what has been modified.
    """

    compressionFilter = tables.Filters(complevel=1, complib='blosc:lz4')
    if pangenome.status["genomesAnnotated"] == "Yes":
        if filename.suffix != ".h5":
            filename = filename.with_suffix(".h5")
        h5f = tables.open_file(filename,"w", filters=compressionFilter)
        logging.getLogger().info("Writing genome annotations...")
        writeAnnotations(pangenome, h5f)
        logging.getLogger().info("Done writing genome annotations")
    elif pangenome.status["genomesAnnotated"] == "Loaded":
        h5f = tables.open_file(filename,"a", filters=compressionFilter)
    else:#if the pangenome is not Yes not Loaded, it's probably not really in a good state ( or something new was coded).
        raise NotImplementedError("Something REALLY unexpected and unplanned for happened here. Dev's contact is ggautrea [at] genoscope [dot] cns [dot] fr.")
    if pangenome.status["genesClustered"] == "Yes":
        logging.getLogger().info("Writing gene families...")
        writeGeneFamilies(pangenome, h5f)
        if pangenome.status["genomesAnnotated"] == "Loaded":
            raise NotImplementedError()
    
    if pangenome.status["NeighborsGraph"] == "Yes":
        writeGraph(pangenome, h5f)
    
    if pangenome.status["Partitionned"] == "Yes":
        raise NotImplementedError()
        ##update geneFamilies with their partition.

    writeStatus(pangenome, h5f)
    h5f.close()
    logging.getLogger().info(f"Done writing the pangenome. It is in file : {filename}")