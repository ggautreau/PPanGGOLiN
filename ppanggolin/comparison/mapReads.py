#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#default libraries
import json
import sys
import os # for size files
import itertools
from multiprocessing import Pool
import glob
import subprocess
from collections import defaultdict

# installed libraries
import numpy as np
import logging # log.txt with times
import tqdm

#local libraries
from ppanggolin.utils import read_compressed_or_not, write_compressed_or_not, get_file_size, filter_fastq, reverse_sign
from ppanggolin.pangenome import Pangenome, GeneFamily, VGPath, Node
from ppanggolin.formats.writeSequences import write_fasta_gene_fam, write_gene_sequences_from_annotations
from ppanggolin.formats.readBinaries import get_gene_sequences_from_file, check_pangenome_info

def map_on_graphs(pangenome, comparison, output_dir, cpu = 1, min_ident = 0.95, disable_bar=False):
    """ Perform a mapping of reads against the pangenome families

    :param pangenome: a pangenome  
    :type pangenome: :class:`ppanggolin.Pangenome`
    :param comparison: a Comparison object filed with samples without counts
    :type comparison: Comparison
    :param output_dir: an output path dir to write mapping results
    :type output_dir: str
    :param cpu: number of cpu to use
    :type cpu: int
    :param min_ident: accept read only if the alignment identity is above this minimum
    :type min_ident: float
    :param disable_bar: either to disable progress bar of not
    :type disable_bar: bool
    :return: a dictionary where families are keys and the values are dicts where samples are keys and values are reads counting
    :rtype: dict
    """
    samples_to_be_mapped = comparison.get_all_samples_having_counts(reversed=True)
    check_pangenome_info(pangenome, need_gene_sequences=True)
    output_local_graphs=output_dir.joinpath("local_graphs/")
    output_combined_graph_and_mapping=output_dir.joinpath("graph_mapping/")
    if not os.path.exists(output_local_graphs):
        os.makedirs(output_local_graphs)
    if not os.path.exists(output_combined_graph_and_mapping):
        os.makedirs(output_combined_graph_and_mapping)

    if not (os.path.exists(output_combined_graph_and_mapping.joinpath("all_graphs.xg")) and
            os.path.exists(output_combined_graph_and_mapping.joinpath("all_graphs.gcsa")) and
            os.path.exists(output_combined_graph_and_mapping.joinpath("all_graphs.snarls")) and
            os.path.exists(output_combined_graph_and_mapping.joinpath("all_graphs.vg")) and
            os.path.exists(output_combined_graph_and_mapping.joinpath("all_graphs.gfa"))):
        make_graphs(pangenome=pangenome, output_local_graphs=output_local_graphs, output_combined_graph_and_mapping=output_combined_graph_and_mapping, cpu=cpu)
    if pangenome.status["variation_graphs"] == "No":
        pangenome.fill_from_GFA_file(output_combined_graph_and_mapping.joinpath("all_graphs.gfa"))

    fasta_nuc_pangenome = output_combined_graph_and_mapping.joinpath("all_nucleotide_families.fasta")
    write_fasta_gene_fam(pangenome, output_combined_graph_and_mapping, "all", 0.95, False, disable_bar=True)

    for sample in tqdm.tqdm(samples_to_be_mapped, disable=disable_bar):
        logging.getLogger().info("filtering by mapping reads over family representatives using minimap2...")
        if sample.is_pair_end():
            #cmd = ["bowtie2", "-p", str(cpu),"--very-sensitive-local", "-x", index_files_prefix, "-1", readFiles[0], "-2", readFiles[1], "--al-conc", output_local_graphs+"/"+IDsample+"_al_conc.fastq", "--al", output_local_graphs+"/"+IDsample+"_al.fastq"]
            cmd = "minimap2 -x sr -t "+ str(cpu)+" "+str(fasta_nuc_pangenome)+" "+str(sample.path_read1)+" "+str(sample.path_read2)+" | cut -f1"
        else:
            cmd = "minimap2 -x sr -t "+ str(cpu)+" "+str(fasta_nuc_pangenome)+" "+str(sample.path_read1) +" | cut -f1"
            #cmd = ["bowtie2", "-p", str(cpu),"--very-sensitive-local", "-x", index_files_prefix, "-U", readFiles[0], "--al", output_local_graphs+"/"+IDsample+"_al.fastq"]
        print(cmd)

        if not (os.path.exists(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq")) or (os.path.exists(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq")) and os.path.exists(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_2.fastq")))):
            selected_reads = str(subprocess.check_output(cmd, shell=True))
            selected_reads = selected_reads[2:]
            selected_reads = selected_reads[:-3]
            selected_reads = selected_reads.split("\\n")
            selected_reads = set(selected_reads)
            filter_fastq(sample.path_read1, output_combined_graph_and_mapping.joinpath(sample.name + "_selected_1.fastq"), selected_reads)
            if sample.is_pair_end():
                filter_fastq(sample.path_read2, output_combined_graph_and_mapping.joinpath(sample.name + "_selected_2.fastq"), selected_reads)
        logging.getLogger().info("mapping reads over the variation graph of each gene family...")
        multimap=False
        if multimap:
            if sample.is_pair_end():
                subprocess.run("vg mpmap -x "+str(output_combined_graph_and_mapping.joinpath("all_graphs.xg"))+ "-g "+str(output_combined_graph_and_mapping.joinpath("all_graphs.gcsa"))+ " -s "+str(output_combined_graph_and_mapping("all_graphs.snarls")) + "-f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq"))+" -f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_2.fastq")) +" -t "+str(cpu)+" -M 10 -m -L 0 | vg view -K -j - > "+str(output_combined_graph_and_mapping.joinpath("mapping.json")), stdout=subprocess.PIPE, shell=True)
            else:
                  subprocess.run("vg mpmap -x "+str(output_combined_graph_and_mapping.joinpath("all_graphs.xg"))+" -g "+str(output_combined_graph_and_mapping.joinpath("all_graphs.gcsa")) +" -s "+str(output_combined_graph_and_mapping("all_graphs.snarls"))+" -f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq"))+" -t "+str(cpu)+" -M 10 -m -L 0 | vg view -K -j - > "+str(output_combined_graph_and_mapping.joinpath("mapping_"+sample.name+".json")), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        else:
            if sample.is_pair_end():
                subprocess.run("vg map --min-ident "+str(min_ident)+" -L 0 -d "+str(output_combined_graph_and_mapping.joinpath("all_graphs"))+" -f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq"))+" -f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_2.fastq"))+" -t "+str(cpu)+" | vg view -a -k - | vg view -K -j - > "+str(output_combined_graph_and_mapping.joinpath("mapping_"+sample.name+".json")), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            else:
                subprocess.run("vg map --min-ident "+str(min_ident)+" -L 0 -d "+str(output_combined_graph_and_mapping.joinpath("all_graphs"))+" -f "+str(output_combined_graph_and_mapping.joinpath(sample.name+"_selected_1.fastq"))+" -t "+str(cpu)+" | vg view -a -k - | vg view -K -j - > "+str(output_combined_graph_and_mapping.joinpath("mapping_"+sample.name+".json")), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        logging.getLogger().info("extracts alignments...")
        pangenome_sample_mapping = pangenome
        pangenome_sample_mapping.__class__ = PangenomeSampleMapping
        pangenome_sample_mapping.extends(sample)
        parse_vgmpmap(output_combined_graph_and_mapping.joinpath("mapping_"+str(sample.name)+".json"), pangenome_sample_mapping, thr=0.7, force_se=False)
        for f in pangenome.gene_families:
            comparison.add_sample_counts_gene_families(sample,f,0)# initialize to zero
        for path in pangenome_sample_mapping.paths:
            count = np.mean([mapped_abundance for mapped_abundance in map(sum, zip(path.unique_mapped_abundances, path.multiple_mapped_abundances))])
            comparison.add_sample_counts_gene_families(sample, path.family, count)

def make_graphs(pangenome, output_local_graphs, output_combined_graph_and_mapping, cpu = 1, disable_bar= False):
    """ Make the variation graphs for all the gene families 

    :param pangenome: a pangenome  
    :type pangenome: :class:`ppanggolin.Pangenome`
    :param output_local_graphs: an output path dir to write locals graphs which correpond to a huge amoumt of dir and files
    :type output_local_graphs: str
    :param output_combined_graph_and_mapping: an output path dir to write combined_graph
    :type output_local_graphs: str
    :param cpu: number of cpu to use
    :type cpu: int
    :param disable_bar: either to disable progress bar of not
    :type disable_bar: bool
    """

    res=[]
    sys.setrecursionlimit(200000)
    args_local_graphs =[]
    
    import copy
    for fam in pangenome.gene_families:
        temp_dir=output_local_graphs.joinpath(str(fam.ID))
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        fam_genes_fasta_path = f"{temp_dir}/{fam.ID}_genes.fasta"
        logging.getLogger().debug(fam_genes_fasta_path)
        with write_compressed_or_not(fam_genes_fasta_path, False) as fasta:
            if pangenome.status["geneSequences"] in ["inFile"]:
                get_gene_sequences_from_file(pangenome.file, fasta, set([gene.ID for gene in fam.genes]),
                                             disable_bar=True)
            elif pangenome.status["geneSequences"] in ["Computed", "Loaded"]:
                write_gene_sequences_from_annotations(fam.genes, fasta, disable_bar=True)
            else:
                # this should never happen if the pangenome has been properly checked before launching this function.
                raise Exception("The pangenome does not include gene sequences")
        args_local_graphs.append([copy.copy(str(fam.ID)), copy.copy(len(fam)), copy.copy(temp_dir), copy.copy(fam_genes_fasta_path)])

    if cpu > 1:
        with Pool(processes = cpu) as p:
            res = [r for r in tqdm.tqdm(p.imap_unordered(make_local_graph, args_local_graphs), disable=disable_bar)]
    else:#for the case where it is called in a daemonic subprocess with a single cpu
        for args_local_graph in tqdm.tqdm(args_local_graphs, disable=disable_bar):
            res.append(make_local_graph(args_local_graph))

    def chunk_list(input_list, chunk_size):
        return [input_list[i:i + chunk_size] for i in range(0, len(input_list), chunk_size)]
    #subprocess.run(['vg', 'ids', '-j', '-c'] + glob.glob(output_local_graphs+"/*/*.vg"))
    nb_chunk=0
    for i, chunk in enumerate(chunk_list(glob.glob(str(output_local_graphs) + "/*/*.vg"), 10000)):
        with open(output_combined_graph_and_mapping.joinpath("all_graphs_chunk_"+str(i)+".vg"),"w") as out_file:
            logging.getLogger().debug(['vg', 'combine'] + glob.glob(str(output_local_graphs)+"/*/*.vg"))
            subprocess.run(['vg', 'combine'] + chunk, stdout=out_file)
        nb_chunk+=1
    with open(output_combined_graph_and_mapping.joinpath("all_graphs.vg"), "w") as out_file:
        logging.getLogger().debug(['vg', 'combine'] + [output_combined_graph_and_mapping.joinpath("all_graphs_chunk_"+str(j)+".vg") for j in range(nb_chunk)])
        subprocess.run(['vg', 'combine'] + [output_combined_graph_and_mapping.joinpath("all_graphs_chunk_"+str(j)+".vg") for j in range(nb_chunk)], stdout=out_file)

    logging.getLogger().debug(['vg', 'index', "-t", str(cpu), "-x", output_combined_graph_and_mapping.joinpath("all_graphs.xg"), output_combined_graph_and_mapping.joinpath("all_graphs.vg")])
    subprocess.run(['vg', 'index', "-t", str(cpu), "-x", output_combined_graph_and_mapping.joinpath("all_graphs.xg"), output_combined_graph_and_mapping.joinpath("all_graphs.vg")])
    logging.getLogger().debug("vg prune -t "+str(cpu)+" "+str(output_combined_graph_and_mapping.joinpath("all_graphs.vg")) +" | vg index -t "+str(cpu)+" -g "+str(output_combined_graph_and_mapping.joinpath("all_graphs.gcsa")) + " - ")
    subprocess.run("vg prune -t "+str(cpu)+" "+str(output_combined_graph_and_mapping.joinpath("all_graphs.vg")) + " | vg index -t "+str(cpu)+" -g "+str(output_combined_graph_and_mapping.joinpath("all_graphs.gcsa")) + " - ", shell=True)
    logging.getLogger().debug("vg snarls -t "+str(cpu)+" "+str(output_combined_graph_and_mapping.joinpath("all_graphs.vg")) + " > "+str(output_combined_graph_and_mapping.joinpath("all_graphs.snarls")))
    subprocess.run("vg snarls -t "+str(cpu)+" "+str(output_combined_graph_and_mapping.joinpath("all_graphs.vg")) + " > " + str(output_combined_graph_and_mapping.joinpath("all_graphs.snarls")), shell=True)
    with open(output_combined_graph_and_mapping.joinpath("all_graphs.gfa"),"w") as out_file:
        logging.getLogger().debug(["vg", "view", "-t", str(cpu), output_combined_graph_and_mapping.joinpath("all_graphs.vg")])
        subprocess.run(["vg", "view", "--threads", str(cpu), output_combined_graph_and_mapping.joinpath("all_graphs.vg")] , stdout=out_file)

def make_local_graph(args):
    fam, nb_genes, temp_dir, fam_genes_fasta_path = args
    logging.getLogger().debug("make_local_graph" + fam)
    cpu = str(1)

    if nb_genes == 1: # if only one sequence in the family, just build a linear graph with vg construct
        logging.getLogger().debug(f"vg construct -t {cpu} -r {fam_genes_fasta_path} -m 256 > {temp_dir}/{fam}.vg")
        subprocess.run(f"vg construct -t {cpu} -r {fam_genes_fasta_path} -m 256 > {temp_dir}/{fam}.vg",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    else:
        logging.getLogger().debug(f"minimap2 -cx asm20 -X -t {cpu} {fam_genes_fasta_path} {fam_genes_fasta_path} > {temp_dir}/{fam}_family_temp.paf")
        subprocess.run(f"minimap2 -cx asm20 -X -t {cpu} {fam_genes_fasta_path} {fam_genes_fasta_path} > {temp_dir}/{fam}_family_temp.paf",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
        logging.getLogger().debug(f"seqwish -t {cpu} -s  {fam_genes_fasta_path} -p {temp_dir}/{fam}_family_temp.paf -b {temp_dir}/ -g {temp_dir}/{fam}_family_temp.gfa")
        subprocess.run(f"seqwish -t {cpu} -s  {fam_genes_fasta_path} -p {temp_dir}/{fam}_family_temp.paf -b {temp_dir}/ -g {temp_dir}/{fam}_family_temp.gfa",shell=True)
        logging.getLogger().debug(f"vg view --threads {cpu} -Fv {temp_dir}/{fam}_family_temp.gfa | vg mod -t {cpu} -n -X 256 -> {temp_dir}/{fam}.vg")
        subprocess.run(f"vg view --threads {cpu} -Fv {temp_dir}/{fam}_family_temp.gfa | vg mod -t {cpu} -n -X 256 -> {temp_dir}/{fam}.vg",shell=True)
        # -n can mess up the graph, check for its integrity otherwise redo the graph without -n
        p1 = subprocess.Popen(["vg","validate",f"{temp_dir}/{fam}.vg"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p1.communicate()
        if err.decode() != "graph: valid\n":
            logging.getLogger().error("error: "+str(err))
            return(False)
    return(True)

class Alignment:

    def __init__(self):
        self.mapped_nodes_id = [] # int ids of nodes
        self.mapped_nodes_cov = [] # coverage of nodes
        self.len = 0
        self.score = 0
        self.identity = 0
        self.nb_aligned = 0
        self.nb_match = 0
        self.nb_subst = 0
        self.nb_indel = 0
        self.nb_errors = 0

def metaalign2align(pangenome: Pangenome, subpath: dict, align_info: dict):
    """
    CONVERT A "META-ALIGNMENT" INTO AN ACTUAL ALIGNMENT PATH
    Alignments in the mapping file (json) are organized into "meta-nodes" (super-nodes containing the actual nodes of the pangenome graph),
    that corresponds to local alignments.
    This function is used as the final step of the BFS function to convert the best (meta)alignments into actual alignments
    (successive nodes id from the pangenome graph) and to compute node-level abundances and errors/variations.
    """
    alignment = Alignment()
    iter_nodes = [n for metanode in align_info["path"] for n in subpath[ metanode ]['path']['mapping']]
    for node in iter_nodes: # for each node get abundance
        nodeID = int(node['position']['node_id'])
        nb_match = 0
        nb_aligned = 0
        for edit in node["edit"]:
            from_len = int(edit.get("from_length", 0))
            to_len = int(edit.get("to_length", 0))
            alignment.nb_indel += abs(from_len-to_len)
            nb_aligned += min(from_len,to_len)
            nb_match += to_len if from_len == to_len and "sequence" not in edit else 0
        alignment.nb_aligned += nb_aligned
        alignment.nb_match += nb_match
        alignment.nb_subst += nb_aligned - nb_match

        abund = nb_aligned/pangenome.nodes[nodeID].len_sequence
        # if node already exists in the list, is the last item and its coverage is not 1 yet, just add the abundance to it (don't duplicate the node)
        if alignment.mapped_nodes_id and alignment.mapped_nodes_id[-1] == nodeID and alignment.mapped_nodes_cov[-1] < 1:
            alignment.mapped_nodes_cov[-1] += abund
        else:
            alignment.mapped_nodes_id.append(nodeID)
            alignment.mapped_nodes_cov.append(abund)
    alignment.len = align_info["read_len"]
    alignment.score = float(align_info["score"]/align_info["read_len"]) # normalized mapping score
    alignment.identity = float(alignment.nb_match/alignment.nb_aligned)
    alignment.nb_errors = float((alignment.nb_indel+alignment.nb_subst)/align_info["read_len"])
    alignment.nb_subst = float(alignment.nb_subst/align_info["read_len"])
    alignment.nb_indel = float(alignment.nb_indel/align_info["read_len"])
    return alignment                

def BFS(aln, pangenome: Pangenome):
    
    subpath = aln['subpath']
    starts = aln["start"] # list of starting meta-nodes (alignment tree nodes that contains pangenome graph nodes)
    
    # step 1: initialize queue
    queue = [*starts]          # contains a set of nodes not yet traversed. 
    
    current_paths = {}  # in runtime, contains all reached paths. key = last path node id, value = tuple(path, read_len, score). 
                        # This enables to quickly retreive all paths that finish on a given node and to conserve the one(s) that minimize the alignment score 
    
    for start in starts: ## there can be several roots, loop over them
        align_info = {}
        align_info["path"] = [start] # list of subpath nodes describing the path of alignement
        align_info["read_len"] = sum([int(edit.get("to_length", 0)) for node in subpath[start]['path']['mapping'] for edit in node["edit"]]) # maybe the alignment can end with a node which has successor, we need to keep track of the position in the read along the alignment
        align_info["score"] = subpath[start].get("score", 0)                                    # score field is not displayed if = 0
        current_paths[start] = [align_info]
    
    ### WALK THE WHOLE SET OF NODES. 
    ### CONSERVE FOR EACH NODE THE BEST PATH (or PATHS in case of equality) SCORE
    while queue:
        
        current_metanode = queue.pop(0) ## pop first element of the queue 
        
        ### END OF THE PROCESS
        if 'next' not in subpath[current_metanode] or all(i["read_len"]==len(aln['sequence']) for i in current_paths[current_metanode]): ## final paths obtained when there is no child or the length of the read is reached
            continue
        
        ### INSPECT ALL CHILDREN
        for child in subpath[current_metanode]['next']:
            
            ### ADD CHILD NODE IF NOT ALREADY IN THE QUEUE: 
            if child not in queue:                                              
                queue.append(child)
            
            child_info = {}
            ### keeping track of the position in the alignment
            child_info["len"] = sum([int(edit.get("to_length", 0)) for node in subpath[child]['path']['mapping'] for edit in node["edit"]]) # maybe the alignment can end with a node which has successor, we need to keep track of the position in the read along the alignment
            ### get child score (or 0 if attribute absent)
            child_info["score"] = subpath[child].get("score", 0)   
                
            ### UPDATE THE SET OF PATHS 
            # 1/ For all detected paths that end with the current_metanode (retreived thanks to `current_paths`)
            # 2/ add all new paths with their scores
            # 3/ remove updated paths from the dictionary
            # 4/ conserve only the best path(s)
            if child not in current_paths:
                current_paths[child]=[]                                                 # init all paths/scores ending with the current node
            for align_info in current_paths[current_metanode]:                              # add all new paths with their scores
                updated_info = {}
                updated_info["path"] = align_info["path"]+[child]
                updated_info["read_len"] = align_info["read_len"] + child_info["len"]
                updated_info["score"] = align_info["score"] + child_info["score"]
                if len(current_paths[child])==0:
                    current_paths[child].append(updated_info)                           # create a new path that finish with the child node
                    continue                                                            # no need to update the paths ending with child, this is the first. 
                ### conserve only the best path(s)
                best_previous_score = current_paths[child][0]["score"]
                if updated_info["score"] > best_previous_score:                         # The new path is better than all previous scores. Create a new set
                    current_paths[child]=[updated_info]                                 # create a new path that finish with the child node
                    continue
                if updated_info["score"] == best_previous_score:                        # The new path is equal to all previous scores, add the new path
                    current_paths[child].append(updated_info)                           # create a new path that finish with the child node
                    continue
                #if updated_info["score"] < best_previous_score: do nothing             # The newly created path is worse than all previous ending on child. Nothing to be done

        current_paths.pop(current_metanode)                                             # No more path end with the current node                                          
    # need an extra-step to check the best scores if paths don't end with the same node
    if len(current_paths) > 1:
        max_score = max([current_paths[end_node][0]['score'] for end_node in current_paths])
        paths_to_remove = [end_node for end_node in current_paths if current_paths[end_node][0]['score'] < max_score]
        [current_paths.pop(end_node) for end_node in paths_to_remove]
    # convert dict to list
    current_paths = list(itertools.chain.from_iterable(current_paths.values()))
    # convert meta-alignment into alignment
    final_paths = list(map(lambda x: metaalign2align(pangenome,subpath,x), current_paths))
    return final_paths

def get_all_alignments_one_read(json_file, pangenome: Pangenome, thr=0.95):
    """ 
    Singled-end: given a first read line, returns all alignemnts correponding to this read.
    Sometimes a read may occur on several successive lines, hence we concatenate the 
    corresponding alignments

    Paired-end: given a first read line, returns all alignemnts correponding to this read and its pair.
    If "name" and "paired_read_name" (mostly cases of distinct read files as opposed to interleaved)are identical, 
    both reads of the pair need to be distinguished based on their sequence.
    """
    line = json_file.readline()

    # End of file EOF
    if not line:
        return None
    
    aln = json.loads(line)
    starting_read_name = aln['name']
    starting_read_pair = aln.get('paired_read_name','')
    mapped_paths = {} # for each read pair, its sequence and list of mapped_paths
    mapped_paths["r1"] = {}
    mapped_paths["r1"]['name'] = aln['name'] # to delete after test phase
    mapped_paths["r1"]["sequence"] = aln['sequence']
    mapped_paths["r1"]["paths"] = []

    # parse only if the read has mapped
    if "subpath" in aln: 
        # get best path(s) for the alignement
        current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
        # store the mapped paths if the score is higher or equal to the threshold
        if  current_mapped_paths[0].score >= thr: 
            mapped_paths["r1"]["paths"] = current_mapped_paths

    while True:

        next_line_to_read_when_leaving_this_function = json_file.tell() 
        line = json_file.readline()
        if not line: 
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break

        aln = json.loads(line)
        if aln['name'] != starting_read_name and aln['name'] != starting_read_pair: # next read
            json_file.seek(next_line_to_read_when_leaving_this_function)
            break

        # parse only if the read has mapped
        if "subpath" in aln: 

            # if first time seeing new sequence, paired-end case
            if aln['sequence'] != mapped_paths["r1"]["sequence"] and "r2" not in mapped_paths:
                mapped_paths["r2"] = {}
                mapped_paths["r2"]['name'] = aln['name'] # to delete after test phase
                mapped_paths["r2"]["sequence"] = aln['sequence']
                mapped_paths["r2"]["paths"] = []

            # get best path(s) for the alignement
            current_mapped_paths = BFS(aln, pangenome) # needs pangenome only for getting node lengths
            # final_paths is a list of Alignments [Alignments]
            current_best_score = current_mapped_paths[0].score

            key = "r1" if aln['sequence'] == mapped_paths["r1"]["sequence"] else "r2"
            if len(mapped_paths[key]["paths"]) == 0: # nothing already higher or equal to thr: 
                if  current_best_score >= thr: 
                    mapped_paths[key]["paths"] += current_mapped_paths
            else: # already something higher or equal to thr: 
                best_stored_score = mapped_paths[key]["paths"][0].score
                if best_stored_score < current_best_score:
                    mapped_paths[key]["paths"] = current_mapped_paths # replace the previously stored paths
                if best_stored_score == current_best_score:
                    mapped_paths[key]["paths"] += current_mapped_paths # add the current paths that have the same score
                # if best_stored_score > current_best_score: do nothing, we do not add those fund paths
    return mapped_paths


class PathSampleMapping(VGPath):
    def extends(self):
        self.unique_mapped_abundances  = [0]*len(self.node_ids) # for each node of the path (ordered as `self.nodes`), store the coverage of mapped reads (each node of a mapped path is set to one, except the two extreme than are usually not 100% covered by the mapped sequence)
        self.total_mapped_unique_reads = 0  #  number of reads with unique mapping on this path.
        self.total_mapped_unique_reads_normalized = 0 # number of coverage ratio reads with unique mapping on this path. Coverage ratio is the length of the read / the len of the sequence of the path
        self.multiple_mapped_abundances = [0]*len(self.node_ids)# for each node of the path (ordered as `self.nodes`), store the coverage of multiple mapped reads (normalized wrt their repartition in other paths)
        self.total_mapped_mult_reads = 0    #  number of reads with corrected multiple mapping on this path.
        self.total_mapped_mult_reads_normalized = 0

class FamilySampleMapping(GeneFamily):
    def extends(self):
        self.unassigned_abundances = {}
        self.reads_name = []  # to delete after test
        self.unassigned_paths = []
        self.unassigned_count = 0
        self.unassigned_count_norm = 0

class PangenomeSampleMapping(Pangenome):
    def extends(self, sample_ID: str):
        self.sample_ID = sample_ID
        self.unassigned_abundances = {}
        self.reads_name = []  # to delete after test
        self.unassigned_paths = []
        self.unassigned_count = 0
        self.unassigned_count_norm = 0

        for f in self.gene_families:
            f.__class__=FamilySampleMapping
            f.extends()

        for p in self.paths:
            p.__class__ = PathSampleMapping
            p.extends()

    def fill_abund_dist(self, gene_match: tuple, aligned_path: Alignment, len_read: int):
        path_id, starting_node_id = gene_match
        nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
        path = self.paths[path_id]
        path.total_mapped_unique_reads += 1
        path.total_mapped_unique_reads_normalized += len_read / self.get_sequence_length(path)
        # Le comptage unique "normalisé" par chemin (incrémentation de (longueur du read)/(longueur du chemin))
        for i in range(nb_mapped_nodes):
            path.unique_mapped_abundances[starting_node_id + i] += aligned_path.mapped_nodes_cov[i]

class NodeSampleMapping(Node):
    def __init__(self, node: Node):
        self.node = node
        self.unique_mapped_abundances = 0

def parse_vgmpmap(json_file_name: str, pangenome_sample_mapping: PangenomeSampleMapping, thr=0.95, force_se=False):
    global NB_READS
    global NB_READS_U
    global NB_READS_M
    global PAIRS
    global SINGLES
    global FILTERED
    global UNASSIGNED
    global ASSIGNED_U
    global ASSIGNED_M
    global INCONSIST
    global INCONSIST_M
    global MULTI_ALIGN
    NB_READS = 0
    NB_READS_U = 0
    NB_READS_M = 0
    PAIRS = 0
    SINGLES = 0
    FILTERED = 0
    UNASSIGNED = 0
    ASSIGNED_U = 0
    ASSIGNED_M = 0
    INCONSIST = 0
    INCONSIST_M = 0
    MULTI_ALIGN = 0

    global SE
    global SAME_CP
    global DIFF_CP
    SE = 0
    SAME_CP = 0
    DIFF_CP = 0

    global ONE_NOALIGN
    global ONE_NOCP
    global AMBIG_CP_RESOLVED
    global AMBIG_RESOLVED
    global AMBIG_REDUCED
    ONE_NOALIGN = 0
    ONE_NOCP = 0
    AMBIG_CP_RESOLVED = 0
    AMBIG_RESOLVED = 0
    AMBIG_REDUCED = 0

    global SE_M
    SE_M = 0

    global SECOND_PASS
    SECOND_PASS = 0

    """
    PARSE MAPPING JSON FILE
    First check all alignments from the same read
    Ignore read if no alignment score > thr
    Ignore read if multiple alignment score > thr

    score = scoring done by vg considering bonus for matches and penalty for mismatches and gap
    identity = Portion of aligned bases that are perfect matches, or 0 if no bases are aligned.
    errors = nb of non aligned bases
    """

    # Optimization: we detect positions in the file of reads with unique mapping. Thus they are not tested twice
    recompute_read = set()

    # DO TWICE THE JOB: Once for detecting the abundance of unique mapped reads
    # FIRST PASS/
    # For each read that maps uniquely: fill the abundance of
    #  1/ each node of the mapped path (extremities are increased <= 1 for each mapped read)
    #  2/ store the abundance of each path simply in term of fully mapped reads
    print("Parsing Alignment: first pass")
    steps = 0
    with open(json_file_name, 'r') as json_file:
        size_file = get_file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell()

            mapped_paths = get_all_alignments_one_read(json_file, pangenome_sample_mapping, thr)

            # end of file
            if mapped_paths == None:
                recompute_read.add(current_seek)  # needed to enable the break during second pass
                break

            NB_READS += len(mapped_paths)
            # paired-end case
            if not force_se and len(mapped_paths) == 2:
                PAIRS += 1

                # 1. no align found for both reads = FILTERED
                if len(mapped_paths["r1"]["paths"]) == 0 and len(mapped_paths["r2"]["paths"]) == 0:
                    FILTERED += 2

                # 1. one of the read has no align = go to single-end style
                elif len(mapped_paths["r1"]["paths"]) == 0 or len(mapped_paths["r2"]["paths"]) == 0:
                    FILTERED += 1
                    ONE_NOALIGN += 1
                    key_to_remove = "r1" if len(mapped_paths["r1"]["paths"]) == 0 else "r2"
                    if key_to_remove == "r1":
                        mapped_paths["r1"] = mapped_paths.pop("r2")
                    else:
                        del mapped_paths["r2"]

                # 1. both reads has alignments
                else:

                    # For each read, parse each alignment path and get corresponding colored path
                    for key in mapped_paths:
                        mapped_paths[key]['colored_match'] = []

                        for i, aligned_path in enumerate(mapped_paths[key]["paths"]):
                            match_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id)
                            mapped_paths[key]["colored_match"] += match_paths
                            for m in match_paths:

                                if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                mapped_paths[key][m[0]].append(i)
                            # get_matching_path for reverse path
                            if len(aligned_path.mapped_nodes_id) > 1:  # don't duplicate the result if the path has only one node
                                match_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id[::-1])
                                mapped_paths[key]["colored_match"] += match_paths
                                for m in match_paths:
                                    if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                    mapped_paths[key][m[0]].append(i)

                    # 2. no cp found for both reads = UNASSIGNED
                    if len(mapped_paths["r1"]["colored_match"]) == 0 and len(mapped_paths["r2"]["colored_match"]) == 0:
                        UNASSIGNED += 2

                        # fill families with unassigned reads (but alignments with several paths possibilities for the same family are not treated)
                        # no shared counts for families as it would add noise for the incompatibility check afterwards
                        for key in mapped_paths:
                            mapped_paths[key]["families"] = [
                                pangenome_sample_mapping.paths[pangenome_sample_mapping.nodes[p.mapped_nodes_id[0]].traversed_path[0]].family.name for
                                p in mapped_paths[key]["paths"]]
                        intersect_families = list(
                            set(mapped_paths["r1"]["families"]) & set(mapped_paths["r2"]["families"]))
                        if len(intersect_families) <= 1:  # if len(intersect_families) > 1, multiple family mapping, do nothing
                            if len(intersect_families) == 1:  # update list_clstr
                                for key in mapped_paths:
                                    mapped_paths[key]["families"] = [x for x in mapped_paths[key]["families"] if
                                                                     x in intersect_families]
                            for key in mapped_paths:
                                if len(mapped_paths[key][
                                           "families"]) == 1:  # if >1, do nothing + if duplicates, do nothing / if ==1, means there is only one path => OK
                                    traversed_family = mapped_paths[key]["families"][0]
                                    pangenome_sample_mapping.get_gene_family(traversed_family).reads_name.append(
                                        mapped_paths[key]["name"])  # to delete after test
                                    pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_paths.append(
                                        mapped_paths[key]["paths"][0].mapped_nodes_id)
                                    pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_count += 1
                                    pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_count_norm += 1 / \
                                                                                                   mapped_paths[key][
                                                                                                       "paths"][0].len
                                    for i, n in enumerate(mapped_paths[key]["paths"][0].mapped_nodes_id):
                                        if n not in pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances:
                                            pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances[n] = 0
                                        pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances[n] += \
                                        mapped_paths[key]["paths"][0].mapped_nodes_cov[i]

                    # 2. one of the read has no cp = go to single-end style
                    # No families intersection here as nothing can be concluded from it
                    elif len(mapped_paths["r1"]["colored_match"]) == 0 or len(mapped_paths["r2"]["colored_match"]) == 0:
                        UNASSIGNED += 1
                        ONE_NOCP += 1
                        key_to_remove = "r1" if len(mapped_paths["r1"]["colored_match"]) == 0 else "r2"

                        # fill families with unassigned reads
                        traversed_families = [
                            pangenome_sample_mapping.paths[pangenome_sample_mapping.nodes[p.mapped_nodes_id[0]].traversed_path[0]].family.name for p in
                            mapped_paths[key_to_remove]["paths"]]
                        if len(traversed_families) == 1:  # if >1, do nothing + if duplicates, do nothing
                            traversed_family = pangenome_sample_mapping.get_gene_family(traversed_families[0])
                            traversed_family.reads_name.append(mapped_paths[key_to_remove]["name"])  # to delete after test
                            traversed_family.unassigned_paths.append(mapped_paths[key_to_remove]["paths"][0].mapped_nodes_id)
                            traversed_family.unassigned_count += 1
                            traversed_family.unassigned_count_norm += 1 / mapped_paths[key_to_remove]["paths"][0].len
                            for i, n in enumerate(mapped_paths[key_to_remove]["paths"][0].mapped_nodes_id):
                                if n not in traversed_family.unassigned_abundances:
                                    traversed_family.unassigned_abundances[n] = 0
                                traversed_family.unassigned_abundances[n] += mapped_paths[key_to_remove]["paths"][0].mapped_nodes_cov[i]

                        # remove key from unassigned read and go to single-end style
                        if key_to_remove == "r1":
                            mapped_paths["r1"] = mapped_paths.pop("r2")
                        else:
                            del mapped_paths["r2"]

                    # 2. both reads has cp
                    else:

                        # colored paths matches intersection
                        list_path_id1 = [cp[0] for cp in mapped_paths['r1']['colored_match']]
                        list_path_id2 = [cp[0] for cp in mapped_paths['r2']['colored_match']]
                        colored_paths_intersect = list(set(list_path_id1) & set(list_path_id2))

                        # 3. length intersection > 1 = multimapping
                        if len(colored_paths_intersect) > 1:
                            recompute_read.add(current_seek)
                            SECOND_PASS += 2
                        # 3. length intersection ==1 = no (more) ambiguity except if duplicates (= multimapping)
                        elif len(colored_paths_intersect) == 1:
                            count1 = list_path_id1.count(colored_paths_intersect[0])
                            count2 = list_path_id2.count(colored_paths_intersect[0])
                            if count1 > 1:
                                MULTI_ALIGN += 1
                            else:
                                idx_cp = [cp[0] for cp in mapped_paths['r1']["colored_match"]].index(
                                    colored_paths_intersect[0])
                                idx_align = mapped_paths['r1'][colored_paths_intersect[0]]
                                ### to delete after test
                                if len(idx_align) > 1:
                                    print("stop")
                                ### end
                                ASSIGNED_U += 1
                                SAME_CP += 1
                                pangenome_sample_mapping.fill_abund_dist(mapped_paths['r1']["colored_match"][idx_cp],
                                                          mapped_paths['r1']["paths"][idx_align[0]],
                                                          len(mapped_paths['r1']["sequence"]))
                                if len(colored_paths_intersect) < len(
                                    mapped_paths['r1']["colored_match"]): AMBIG_CP_RESOLVED += 1
                            if count2 > 1:
                                MULTI_ALIGN += 1
                            else:
                                idx_cp = [cp[0] for cp in mapped_paths['r2']["colored_match"]].index(
                                    colored_paths_intersect[0])
                                idx_align = mapped_paths['r2'][colored_paths_intersect[0]]
                                ### to delete after test
                                if len(idx_align) > 1:
                                    print("stop")
                                ### end
                                ASSIGNED_U += 1
                                SAME_CP += 1
                                pangenome_sample_mapping.fill_abund_dist(mapped_paths['r2']["colored_match"][idx_cp],
                                                          mapped_paths['r2']["paths"][idx_align[0]],
                                                          len(mapped_paths['r2']["sequence"]))
                                if len(colored_paths_intersect) < len(
                                    mapped_paths['r2']["colored_match"]): AMBIG_CP_RESOLVED += 1
                        # 3. no intersection = check strains intersection
                        else:

                            # strains intersection
                            for key in mapped_paths:
                                mapped_paths[key]['corresp_orgs'] = []
                                for i, c in enumerate(mapped_paths[key]["colored_match"]):
                                    corresp_orgs = list(pangenome_sample_mapping.paths[c[0]].organisms.keys())
                                    mapped_paths[key]["corresp_orgs"] += corresp_orgs
                                    for s in corresp_orgs:
                                        if s not in mapped_paths[key]: mapped_paths[key][s] = []
                                        mapped_paths[key][s].append(i)
                            org_intersect = list(
                                set(mapped_paths['r1']['corresp_orgs']) & set(mapped_paths['r2']['corresp_orgs']))

                            # reduce original list of colored path
                            if len(org_intersect) >= 1:

                                for key in mapped_paths:
                                    update_cp = []
                                    for s in org_intersect:
                                        update_cp += [mapped_paths[key]["colored_match"][i] for i in
                                                      mapped_paths[key][s]]
                                    if len(update_cp) < len(mapped_paths[key]["colored_match"]):
                                        if len(update_cp) == 1:
                                            AMBIG_RESOLVED += 1
                                        else:
                                            AMBIG_REDUCED += 1
                                    mapped_paths[key]["colored_match"] = update_cp

                                    if len(mapped_paths[key]["colored_match"]) > 1:
                                        recompute_read.add(current_seek)
                                        SECOND_PASS += 1
                                    else:  # len(mapped_paths[key]["colored_match"]) == 1
                                        corresp_align = [mapped_paths[key]['paths'][i] for i in
                                                         mapped_paths[key][mapped_paths[key]["colored_match"][0][0]]]
                                        if len(corresp_align) == 1:
                                            ASSIGNED_U += 1
                                            DIFF_CP += 1
                                            pangenome_sample_mapping.fill_abund_dist(mapped_paths[key]["colored_match"][0],
                                                                      corresp_align[0],
                                                                      len(mapped_paths[key]["sequence"]))
                                        else:
                                            MULTI_ALIGN += 1
                            else:
                                INCONSIST += 2

            # single-end case
            if len(mapped_paths) == 1 or force_se:

                for key in mapped_paths:

                    SINGLES += 1

                    name = mapped_paths[key]['name']  # to delete after test phase
                    aligned_read = mapped_paths[key]["sequence"]
                    mapped_paths_temp = mapped_paths[key]["paths"]

                    if len(mapped_paths_temp) == 0:
                        FILTERED += 1
                        continue  # no path found

                    if len(mapped_paths_temp) > 1:
                        recompute_read.add(current_seek)
                        SECOND_PASS += 1
                        continue  # Here we deal only with reads mapping exactly one path

                    # if len(mapped_paths_temp) == 1
                    aligned_path = mapped_paths_temp[0]  # for clarity
                    match_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id)
                    if len(aligned_path.mapped_nodes_id) > 1:  # don't duplicate the result if the path has only one node
                        match_paths += pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id[::-1])

                    if len(match_paths) > 1:  # we may have several paths corresponding to a unique alignment
                        recompute_read.add(current_seek)
                        SECOND_PASS += 1
                        continue
                    if len(match_paths) == 0:
                        UNASSIGNED += 1
                        # add info to corresponding family
                        traversed_family = pangenome_sample_mapping.paths[pangenome_sample_mapping.nodes[aligned_path.mapped_nodes_id[0]].traversed_path[0]].family.name
                        pangenome_sample_mapping.get_gene_family(traversed_family).reads_name.append(name)  # to delete after test
                        pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_paths.append(aligned_path.mapped_nodes_id)
                        pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_count += 1
                        pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_count_norm += 1 / aligned_path.len
                        for i, n in enumerate(aligned_path.mapped_nodes_id):
                            if n not in pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances:
                                pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances[n] = 0
                            pangenome_sample_mapping.get_gene_family(traversed_family).unassigned_abundances[n] += aligned_path.mapped_nodes_cov[i]
                        continue
                    else:
                        ASSIGNED_U += 1
                        SE += 1
                        continue

                    # if len(match_paths) == 1
                    pangenome_sample_mapping.fill_abund_dist(match_paths[0], aligned_path, len(aligned_read))

    # DO TWICE THE JOB: Once for detecting the abundance of unique mapped reads
    # Once for dealing with multimapped reads
    # SECOND PASS/
    # For each read that maps on several paths
    # detect the total_mapped_unique_reads of each of the mapped paths
    # This provides an abundance a,b,c eg for 3 mapped paths respectively A, B, C.
    # For path 'A', add in each node A.multiple_mapped_abundances[node] a/(a+b+c)
    # For path 'B', add in each node B.multiple_mapped_abundances[node] b/(a+b+c)
    # For path 'C', add in each node C.multiple_mapped_abundances[node] c/(a+b+c)

    steps = 0
    print("Parsing Alignment: second pass")
    with open(json_file_name, 'r') as json_file:
        size_file = get_file_size(json_file)
        while True:
            steps += 1
            current_seek = json_file.tell()
            if current_seek not in recompute_read:
                json_file.readline()  # dont care
                continue
            mapped_paths = get_all_alignments_one_read(json_file, pangenome_sample_mapping, thr)

            # end of file
            if mapped_paths == None:
                break

            # paired-end case
            if not force_se and len(mapped_paths) == 2:

                # 1. one of the read has no align = go to single-end style
                if len(mapped_paths["r1"]["paths"]) == 0 or len(mapped_paths["r2"]["paths"]) == 0:
                    key_to_remove = "r1" if len(mapped_paths["r1"]["paths"]) == 0 else "r2"
                    if key_to_remove == "r1":
                        mapped_paths["r1"] = mapped_paths.pop("r2")
                    else:
                        del mapped_paths["r2"]

                # 1. both reads has alignments
                else:

                    # For each read, parse each alignment path and get corresponding colored path
                    for key in mapped_paths:
                        mapped_paths[key]['colored_match'] = []

                        for i, aligned_path in enumerate(mapped_paths[key]["paths"]):
                            match_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id)
                            mapped_paths[key]["colored_match"] += match_paths
                            for m in match_paths:
                                if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                mapped_paths[key][m[0]].append(i)
                            # get_matching_path for reverse path
                            if len(aligned_path.mapped_nodes_id) > 1:  # don't duplicate the result if the path has only one node
                                match_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id[::-1])
                                mapped_paths[key]["colored_match"] += match_paths
                                for m in match_paths:
                                    if m[0] not in mapped_paths[key]: mapped_paths[key][m[0]] = []
                                    mapped_paths[key][m[0]].append(i)

                    # 2. one of the read has no cp = go to single-end style
                    if len(mapped_paths["r1"]["colored_match"]) == 0 or len(mapped_paths["r2"]["colored_match"]) == 0:
                        key_to_remove = "r1" if len(mapped_paths["r1"]["colored_match"]) == 0 else "r2"
                        if key_to_remove == "r1":
                            mapped_paths["r1"] = mapped_paths.pop("r2")
                        else:
                            del mapped_paths["r2"]
                    # 2. both reads has cp
                    else:

                        # colored paths matches intersection
                        list_path_id1 = [cp[0] for cp in mapped_paths['r1']['colored_match']]
                        list_path_id2 = [cp[0] for cp in mapped_paths['r2']['colored_match']]
                        colored_paths_intersect = list(set(list_path_id1) & set(list_path_id2))

                        # 3. length intersection >= 1 ( ==1 should not be happening in second pass)
                        if len(colored_paths_intersect) >= 1:

                            for key in mapped_paths:

                                # if all final cp has only one align -> OK, else too complex -> MULTI_ALIGN count
                                if all([len(mapped_paths[key][cp]) == 1 for cp in colored_paths_intersect]):

                                    ASSIGNED_M += 1
                                    # compute a+b+c (cf earlier comments)
                                    sum_covered_paths = sum([pangenome_sample_mapping.paths[cp].total_mapped_unique_reads for cp in
                                                             colored_paths_intersect])
                                    # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                                    for cp in mapped_paths[key]['colored_match']:
                                        if cp[0] in colored_paths_intersect:
                                            path_id = cp[0]
                                            starting_node_id = cp[1]
                                            idx_align = mapped_paths[key][cp[0]][0]
                                            aligned_path = mapped_paths[key]['paths'][idx_align]
                                            nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                                            path = pangenome_sample_mapping.paths[path_id]
                                            if sum_covered_paths == 0:  # if no unique mapped reads, equal repartition to the strains
                                                ratio = 1 / len(colored_paths_intersect)
                                            else:
                                                ratio = (path.total_mapped_unique_reads) / float(sum_covered_paths)
                                            path.total_mapped_mult_reads += ratio
                                            path.total_mapped_mult_reads_normalized += ratio * len(
                                                mapped_paths[key]['sequence']) / pangenome_sample_mapping.get_sequence_length(path)
                                            for i in range(nb_mapped_nodes):
                                                path.multiple_mapped_abundances[starting_node_id + i] += \
                                                aligned_path.mapped_nodes_cov[i] * ratio
                                else:
                                    MULTI_ALIGN += 1

                        # 3. no intersection = check strains intersection
                        else:

                            # strains intersection
                            for key in mapped_paths:
                                mapped_paths[key]['corresp_orgs'] = []
                                for i, c in enumerate(mapped_paths[key]["colored_match"]):
                                    corresp_orgs = list(pangenome_sample_mapping.paths[c[0]].organisms.keys())
                                    mapped_paths[key]["corresp_orgs"] += corresp_orgs
                                    for s in corresp_orgs:
                                        if s not in mapped_paths[key]: mapped_paths[key][s] = []
                                        mapped_paths[key][s].append(i)
                            org_intersect = list(
                                set(mapped_paths['r1']['corresp_orgs']) & set(mapped_paths['r2']['corresp_orgs']))

                            # reduce original list of colored path
                            if len(org_intersect) >= 1:

                                for key in mapped_paths:
                                    update_cp = []
                                    for s in org_intersect:
                                        update_cp += [mapped_paths[key]["colored_match"][i] for i in
                                                      mapped_paths[key][s]]
                                    mapped_paths[key]["colored_match"] = update_cp

                                    # if len(mapped_paths[key]["colored_match"]) == 1, read has already been processed in first pass, don't process it twice
                                    if len(mapped_paths[key]["colored_match"]) > 1:

                                        # if all final cp has only one align -> OK, else too complex -> MULTI_ALIGN count
                                        if all([len(mapped_paths[key][cp[0]]) == 1 for cp in
                                                mapped_paths[key]["colored_match"]]):

                                            ASSIGNED_M += 1
                                            # compute a+b+c (cf earlier comments)
                                            sum_covered_paths = sum(
                                                [pangenome_sample_mapping.paths[cp[0]].total_mapped_unique_reads for cp in
                                                 mapped_paths[key]["colored_match"]])
                                            # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                                            for cp in mapped_paths[key]['colored_match']:
                                                path_id = cp[0]
                                                starting_node_id = cp[1]
                                                idx_align = mapped_paths[key][cp[0]][0]
                                                aligned_path = mapped_paths[key]['paths'][idx_align]
                                                nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                                                path = pangenome_sample_mapping.paths[path_id]
                                                if sum_covered_paths == 0:  # if no unique mapped reads, equal repartition to the strains
                                                    ratio = 1 / len(mapped_paths[key]["colored_match"])
                                                else:
                                                    ratio = (path.total_mapped_unique_reads) / float(sum_covered_paths)
                                                path.total_mapped_mult_reads += ratio
                                                path.total_mapped_mult_reads_normalized += ratio * len(
                                                    mapped_paths[key]['sequence']) / pangenome_sample_mapping.get_sequence_length(path)
                                                for i in range(nb_mapped_nodes):
                                                    path.multiple_mapped_abundances[starting_node_id + i] += \
                                                    aligned_path.mapped_nodes_cov[i] * ratio
                                        else:
                                            MULTI_ALIGN += 1
                            else:
                                if len(mapped_paths['r1']["colored_match"]) == 1 or len(
                                        mapped_paths['r2']["colored_match"]) == 1:
                                    INCONSIST_M += 1
                                else:
                                    INCONSIST_M += 2

            # single-end case
            if len(mapped_paths) == 1 or force_se:

                for key in mapped_paths:

                    name = mapped_paths[key]["name"]  # to delete after test phase
                    aligned_read = mapped_paths[key]["sequence"]
                    mapped_paths_temp = mapped_paths[key]["paths"]

                    # we retreive the paths corresponding to this alignments:
                    flag = False
                    for aligned_path in mapped_paths_temp:
                        found_gene_paths = pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id)
                        if len(aligned_path.mapped_nodes_id) > 1:  # don't duplicate the result if the path has only one node
                            found_gene_paths += pangenome_sample_mapping.get_matching_path(aligned_path.mapped_nodes_id[::-1])

                        if len(found_gene_paths) > 0:
                            flag = True

                        # compute a+b+c (cf earlier comments)
                        sum_covered_paths = 0
                        for found_gene_path in found_gene_paths:
                            path_id = found_gene_path[0]
                            path = pangenome_sample_mapping.paths[path_id]
                            sum_covered_paths += path.total_mapped_unique_reads  # TODO: valider avec Kevin ce +1 (en cas de tout à zero)

                        # fill corresponding nodes normalized abundances (a/(a+b+c) cf earlier comments
                        for found_gene_path in found_gene_paths:
                            path_id = found_gene_path[0]
                            starting_node_id = found_gene_path[1]
                            nb_mapped_nodes = len(aligned_path.mapped_nodes_id)
                            path = pangenome_sample_mapping.paths[path_id]
                            if sum_covered_paths == 0:  # if no unique mapped reads, equal repartition to the strains
                                ratio = 1 / len(found_gene_paths)
                            else:
                                ratio = (path.total_mapped_unique_reads) / float(sum_covered_paths)
                            path.total_mapped_mult_reads += ratio
                            path.total_mapped_mult_reads_normalized += ratio * len(
                                aligned_read) / pangenome_sample_mapping.get_sequence_length(path)
                            for i in range(nb_mapped_nodes):
                                path.multiple_mapped_abundances[starting_node_id + i] += aligned_path.mapped_nodes_cov[
                                                                                             i] * ratio

                    if flag:
                        ASSIGNED_M += 1
                        SE_M += 1
                    else:
                        UNASSIGNED += 1
