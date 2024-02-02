#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import pdb
import tempfile
import subprocess
import argparse
from collections import defaultdict
import sys
import os
from multiprocessing import Pool
import textwrap
from operator import xor
from math import inf

# installed libraries
import numpy as np
from scipy.stats import ttest_rel, ttest_ind, wilcoxon, ranksums, iqr, fisher_exact, chi2_contingency, hmean
from scipy.stats.contingency import association
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import StratifiedTable
from scipy.spatial.distance import pdist
from scipy.sparse import csc_matrix
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.sparse import csc_matrix
import plotly.graph_objs as go
import plotly.offline as out_plotly
import plotly.figure_factory as ff
import colorlover as cl

#local libraries
from ppanggolin.utils import mk_outdir, restricted_float# read_compressed_or_not, write_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.genome import Organism
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.comparison.mapReads import map_on_graphs

class Sample:
    def __init__(self, name, path_read1 = None, path_read2 = None):
        self.name = name
        self.path_read1 = path_read1
        self.path_read2 = path_read2
        self.imported_counts = False
        self.sample_family_persistent_coverage = []
        self.covered = False
        self.th = float("nan")

    @property
    def sample_family_persistent_covered(self):
        return(len([c for c in self.sample_family_persistent_coverage if c > 0]))

    def is_pair_end(self):
        if self.path_read1 is not None and self.path_read2 is not None:
            return(True)
        else:
            return(False)

class FamilyComparisonResult:
    def __init__(self, fam):
        self.family = fam
        self.contingency_table = None
        self.pvalue = float("nan")
        self.oddsratio =  float("nan")
        self.cramer_V = float("nan")
        self.corrected_pvalue = float("nan")
        self.module_results = None

class ModuleComparisonResult:
    def __init__(self, module, family_results):
        self.module = module
        self.family_results = family_results

    @property
    def combined_pvalue(self):
        return(hmean([fr.pvalue for fr in self.family_results]))

    @property
    def combined_corrected_pvalue(self):
        return(hmean([fr.corrected_pvalue for fr in self.family_results]))

    @property
    def mean_cramer_V(self):
        return(np.mean([fr.cramer_V for fr in self.family_results]))
    @property
    def mean_oddsratio(self):
        return (np.mean([fr.oddsratio for fr in self.family_results]))
    @property
    def pooled_oddsratio(self):
        return(StratifiedTable([fr.contingency_table for fr in self.family_results]).oddsratio_pooled)

class Comparison:
    def __init__(self, pangenome, condition1_name, condition2_name, name = None):
        self.pangenome = pangenome
        self.condition1_name = condition1_name
        self.condition2_name = condition2_name
        self.condition1_samples_or_orgs = {}
        self.condition2_samples_or_orgs = {}

        self.sample_counts_gene_families = dict(zip(self.pangenome.gene_families,[dict() for i in range(self.pangenome.number_of_gene_families())]))
        self.performed = False
        self.results = dict(zip(self.pangenome.gene_families,[None]*self.pangenome.number_of_gene_families()))
        if name is None:
            self.name = self.condition1_name+"_vs_"+self.condition2_name

    def add_sample_or_org(self, sample_or_org, condition = None):
        if condition is not None:
            if str(condition) == "1":
                self.condition1_samples_or_orgs[sample_or_org.name] = sample_or_org
            if str(condition) == "2":
                self.condition2_samples_or_orgs[sample_or_org.name] = sample_or_org
            self.performed = False
        else:
            raise Exception("No condition provided")

    def get_sample_or_org(self, sample_or_org_name : str):
        if sample_or_org_name in self.condition1_samples_or_orgs:
            return(self.condition1_samples_or_orgs[sample_or_org_name])
        elif sample_or_org_name in self.condition2_samples_or_orgs:
            return(self.condition2_samples_or_orgs[sample_or_org_name])
        else:
            raise KeyError("{sample_or_org_name} is absent of this comparison")

    def add_sample_counts_gene_families(self, sample : Sample, family : GeneFamily, count : float):
        if sample in self.sample_counts_gene_families[family]:
            self.sample_counts_gene_families[family][sample] += count
        else:
            self.sample_counts_gene_families[family][sample] = count
        sample.imported_counts=True
        self.performed = False

    def import_count_matrix(self, count_matrix_file, sep = "\t"):
        """ Import a count matrix corresponding to a previously computed mapping of reads against the pangenome families
        :param count_matrix_file: input path
        :type count_matrix_file: str
        :return: a dictionary where families are keys and the values are dicts where samples are keys and values are reads counting
        :rtype: dict
        """
        sample_names=None
        with open(count_matrix_file,"r") as count_matrix:
            for line in count_matrix:
                elements = line.strip("\n").split(sep)
                if sample_names is None:# header
                    sample_names = elements[1:]
                else:
                    for index_sample, count in enumerate(elements[1:]):
                        self.add_sample_counts_gene_families(self.get_sample_or_org(sample_names[index_sample]),self.pangenome.get_gene_family(elements[0]),float(count))
        self.performed = False

    def get_all_orgs(self, condition = "both"):
        ret = set()
        condition = str(condition)
        if condition == "both" or condition[-1] == "1":
            for samples_or_orgs in self.condition1_samples_or_orgs.values():
                if samples_or_orgs.__class__ == Organism:
                    ret.add(samples_or_orgs)
        if condition == "both" or condition[-1] == "2":
            for samples_or_orgs in self.condition2_samples_or_orgs.values():
                if samples_or_orgs.__class__ == Organism:
                    ret.add(samples_or_orgs)
        return(ret)

    def get_all_samples(self, condition = "both"):
        ret = set()
        condition = str(condition)
        if condition == "both" or condition[-1] == "1":
            for samples_or_orgs in self.condition1_samples_or_orgs.values():
                if samples_or_orgs.__class__ == Sample:
                    ret.add(samples_or_orgs)
        if condition == "both" or condition[-1] == "2":
            for samples_or_orgs in self.condition2_samples_or_orgs.values():
                if samples_or_orgs.__class__ == Sample:
                    ret.add(samples_or_orgs)
        return(ret)

    def get_all_covered_samples(self, reversed=False, condition = "both"):
        ret = set()
        condition = str(condition)
        for s in self.get_all_samples(condition = condition):
            if xor(s.covered, reversed):
                ret.add(s)
        return(ret)

    def get_all_samples_or_orgs(self, condition = "both"):
        ret = set()
        condition = str(condition)
        if condition == "both" or condition[-1] == "1":
            for samples_or_orgs in self.condition1_samples_or_orgs.values():
                ret.add(samples_or_orgs)
        if condition == "both" or condition[-1] == "2":
            for samples_or_orgs in self.condition2_samples_or_orgs.values():
                ret.add(samples_or_orgs)
        return(ret)

    def get_covered_sample_family_presence(self, f):
        return(self.get_all_covered_samples().intersection(set([s for s, count in self.sample_counts_gene_families[f].items() if count > s.th])))

    def get_all_samples_having_counts(self, reversed=False, condition = "both"):
        ret = set()
        for s in self.get_all_samples(condition = condition):
            if xor(s.imported_counts, reversed):
                ret.add(s)
        return(ret)

    def write_conditions(self, dir = None):
        with open(dir + "/conditions_" + self.name + ".tsv", "w") as conditions:
            samples = list(self.get_all_samples(condition=1) | self.get_all_samples(condition=2))
            orgs = list(self.get_all_orgs(condition=1) | self.get_all_orgs(condition=2))
            nb_persistent = 0
            for f in self.pangenome.gene_families:
                if f.partition == "P":
                    nb_persistent += 1
            conditions.write("sample_or_orgs\tcondition\tpersistent_covered\tth_sample\tnb_persistent_covered\tratio_persitent_covered\n")
            for sample_or_org in samples+orgs:
                if sample_or_org.name in self.condition1_samples_or_orgs and sample_or_org.name not in self.condition2_samples_or_orgs:
                    if sample_or_org in self.get_all_samples():
                        nb_covered = sample_or_org.sample_family_persistent_covered
                        conditions.write(f'{sample_or_org.name}\t{self.condition1_name}\t{sample_or_org.covered}\t{sample_or_org.th}\t{nb_covered}\t{nb_covered/nb_persistent}\n')
                    else:
                        conditions.write(sample_or_org.name + "\t"+self.condition1_name+"\t\t\t\t\n")
                elif sample_or_org.name in self.condition2_samples_or_orgs and sample_or_org.name not in self.condition1_samples_or_orgs:
                    if sample_or_org in self.get_all_samples():
                        nb_covered = sample_or_org.sample_family_persistent_covered
                        conditions.write(f'{sample_or_org.name}\t{self.condition2_name}\t{sample_or_org.covered}\t{sample_or_org.th}\t{nb_covered}\t{nb_covered/nb_persistent}\n')
                    else:
                        conditions.write(sample_or_org.name + "\t" + self.condition2_name + "\t\t\t\t\n")

    def write_used_dataset_binary_matrix(self, condition = "both", dir = None):
        with open(dir + "/dataset_" + condition + "_" + self.name + "_binary.tsv", "w") as dataset:
            samples = list(self.get_all_covered_samples(condition=condition))
            orgs = list(self.get_all_orgs(condition=condition))
            dataset.write("\t".join(["gene_families"] + [str(s_or_o.name) for s_or_o in samples + orgs]) + "\n")
            for gf in self.pangenome.gene_families:
                dataset.write("\t".join([str(gf.name)] + ["1" if self.sample_counts_gene_families[gf][s] > s.th else "0" for s in samples] + ["1" if org in gf.organisms else "0" for org in orgs]) + "\n")

    def write_used_dataset_count(self, condition = "both", dir = None):
        with open(dir + "/dataset_" + condition + "_" + self.name + "_counts.tsv", "w") as dataset:
            samples = list(self.get_all_samples(condition=condition).intersection(self.get_all_covered_samples()))
            dataset.write("\t".join(["#th_samples"] + [str(s.th) for s in samples]) + "\n")
            dataset.write("\t".join(["gene_families"] + [str(s.name) for s in samples]) + "\n")
            for gf in self.pangenome.gene_families:
                dataset.write("\t".join(
                    [str(gf.name)] + [str(self.sample_counts_gene_families[gf][s]) for s in samples]) + "\n")

    def perform_comparisons(self, functional_modules = None, min_cov_persistent=0.85, th_ratio_persistent_mean_coverage=0.05, dir = None):
        """
        Perform comparations based on the occurences of genes in a pangenome between two list of genomes corresponding to 2 conditions.
        Reads a pangenome object and return a dictionnary where gene families are the keys and the values are lists of 4 elements (p-value, oddsratio, V-cramer,p-value-corrected).

        :param functional_modules: a dictionary where families are keys and sets of modules are values
        :type functional_modules: dict
        :param min_cov_persistent: The name of different strains for condition1
        :type min_cov_persistent: float
        :param th_ratio_persistent_mean_coverage: The name of different strains for condition2
        :type th_ratio_persistent_mean_coverage: float
        :param dir: output directory path
        :type dir: str
        """
        uncorrected_pvalues = []
        nb_persistent = 0
        for f in self.pangenome.gene_families:
            if f.partition == "P":
                nb_persistent += 1
                for s, abundance in self.sample_counts_gene_families[f].items():
                    if abundance > 0:
                        s.sample_family_persistent_coverage.append(abundance)
        for s in self.get_all_samples():
            if s.sample_family_persistent_covered > nb_persistent * min_cov_persistent:
                s.covered = True
                s.th = sum(s.sample_family_persistent_coverage) / s.sample_family_persistent_covered * th_ratio_persistent_mean_coverage

        for f in self.pangenome.gene_families:
            dataset1_fampresence = 0
            dataset2_fampresence = 0
            dataset1_famabsence = 0
            dataset2_famabsence = 0
            self.results[f] = FamilyComparisonResult(f)
            if len(self.get_all_samples())>0:
                family_in_org_samples = set([org for org in f.organisms]).union(self.get_covered_sample_family_presence(f))
            else:
                family_in_org_samples = set([org for org in f.organisms])
            for org_sample in self.condition1_samples_or_orgs.values():
                if org_sample in family_in_org_samples:
                    dataset1_fampresence += 1
                else:
                    dataset1_famabsence += 1
            for org_sample in self.condition2_samples_or_orgs.values():
                if org_sample in family_in_org_samples:
                    dataset2_fampresence += 1
                else:
                    dataset2_famabsence += 1

            self.results[f].contingency_table = np.array([[dataset1_fampresence, dataset2_fampresence], [dataset1_famabsence, dataset2_famabsence]])
            try:
                self.results[f].oddsratio, self.results[f].pvalue = fisher_exact(self.results[f].contingency_table)
                self.results[f].cramer_V = association(self.results[f].contingency_table, method="cramer")
            except:
                pass

            uncorrected_pvalues.append(self.results[f].pvalue)

        self.write_used_dataset_binary_matrix("both", dir)
        self.write_conditions(dir)
        #cmd = "pyseer --no-distances --phenotypes "+dir+"/conditions_" + self.name + ".tsv --pres "+dir+"dataset_both_" + self.name + "_binary.tsv > "+dir+"/pyseer.assoc"
        #logging.getLogger().info(cmd)
        #subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if dir is not None:
            self.write_used_dataset_binary_matrix("1", dir)
            self.write_used_dataset_binary_matrix("2", dir)
            self.write_used_dataset_count("1", dir)
            self.write_used_dataset_count("2", dir)
        all_corrected_pvals = multipletests(uncorrected_pvalues, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        for index, f in enumerate(self.pangenome.gene_families):
            self.results[f].corrected_pvalue = all_corrected_pvals[index]
        if functional_modules is not None:
            for module, fams in functional_modules.items():
                mr = ModuleComparisonResult(module, [self.results[self.pangenome.get_gene_family(f)] for f in fams])
                for f in fams:
                    self.results[self.pangenome.get_gene_family(f)].module_results = mr
        with open(dir + "/results_" + self.name + ".tsv", "w") as tsvfile:
            tsvfile.write("Id\tGene_family_name\tpartition\tpvalue\toddsratio\tV\tcorrected_pvalue\tmodule_id\tmodule_combined_pvalue\tmodule_pooled_oddsratio\tmodule_mean_cramer_V\tmodule_combined_corrected_pvalue\n")
            for fam_result in sorted(self.results.values(), key=lambda x: x.pvalue):
                index = str(fam_result.family.ID)
                fam_name = fam_result.family.name
                partition = fam_result.family.named_partition
                pvalue_str = str(fam_result.pvalue)
                oddsratio_str = str(fam_result.oddsratio)
                cramer_V_str = str(fam_result.cramer_V)
                corrected_pvalue_str = str(fam_result.corrected_pvalue)
                try:
                    module_id = str(fam_result.module_results.module)
                    module_combined_pvalue_str = str(fam_result.module_results.combined_pvalue)
                    module_pooled_oddsratio_str = str(fam_result.module_results.pooled_oddsratio)
                    module_mean_cramer_V_str = str(fam_result.module_results.mean_cramer_V)
                    module_combined_corrected_pvalue_str = str(fam_result.module_results.combined_corrected_pvalue)
                    tsvfile.write(
                        index + "\t" + fam_name + "\t" + partition + "\t" + pvalue_str + "\t" + oddsratio_str + "\t" + cramer_V_str + "\t" + corrected_pvalue_str + "\t" + module_id + "\t" + module_combined_pvalue_str + "\t" + module_pooled_oddsratio_str + "\t" + module_mean_cramer_V_str + "\t" + module_combined_corrected_pvalue_str + "\n")
                except:
                    tsvfile.write(
                        index + "\t" + fam_name + "\t" + partition + "\t" + pvalue_str + "\t" + oddsratio_str + "\t" + cramer_V_str + "\t" + corrected_pvalue_str + "\t\t\t\t\t\n")

        self.performed = True

    def draw_tileplot_samples(self, output, nocloud=False, filter_by_stats=True):
        sys.setrecursionlimit(100000)
        if len(self.get_all_covered_samples()) > 500 and nocloud is False:
            logging.getLogger().warning("You asked to draw a tile plot for a lot of samples (>500). Your browser will probably not be able to open it.")
        if len(self.get_all_covered_samples()) <2 :
            logging.getLogger().warning("Plotting aborted. Not enough covered samples.")
            return(None)
        # samples_dataset1 = []
        # samples_dataset2 = []
        # for name_c, c in pangenome.comparisons.items():
        #     samples_dataset1 += [s.ID for s in list(c.dataset1.samples_dataset)]
        #     samples_dataset2 += [s.ID for s in list(c.dataset2.samples_dataset)]
        # samples = samples_dataset1 + samples_dataset2

        logging.getLogger().info("Drawing the sample tile plot...")
        data = []
        if nocloud:
            families = {fam for fam in self.pangenome.gene_families if not fam.partition.startswith("C")}
        else:
            families = set(self.pangenome.gene_families)
        sample_index = {v: k for k, v in enumerate(self.get_all_covered_samples())}
        index2sample = {k: v for k, v in enumerate(self.get_all_covered_samples())}
        COLORS = {"pangenome": "black", "exact_accessory": "#EB37ED", "exact_core": "#FF2828", "soft_core": "#c7c938",
                  "soft_accessory": "#996633", "shell": "#00D860", "persistent": "#F7A507", "cloud": "#79DEFF",
                  "undefined": "#828282"}
        logging.getLogger().info("start with matrice")
        persistent_fam = [fam for fam in families if fam.named_partition == "persistent"]
        shell_fam = [fam for fam in families if fam.named_partition == "shell"]
        cloud_fam = [fam for fam in families if fam.named_partition == "cloud"]
        COLORS_partitions = {"shell": "#00D860", "persistent": "#F7A507", "cloud": "#79DEFF", "undefined": "#828282"}
        # for sample in sample_index:
        #     counts_persistent = [comparison.samples[sample].gene_families_map_count[fam.name] for fam in persistent_fam]
        #     counts_shell = [comparison.samples[sample].gene_families_map_count[fam.name] for fam in shell_fam]
        #     # counts_cloud = [pangenome.samples[sample].gene_families_map_count[fam.name] for fam in cloud_fam]
        #     pangenome.samples[sample].th = np.quantile(counts_persistent, 0.01)
        #     # distinguer persistent, shell, cloud, 0
        #     # data.append([1 if pangenome.samples[sample].gene_families_map_count[fam.name]>sample_th[sample] else 0 for fam in families])
        #     import os
        #     import plotly.figure_factory as ff
        #     if not os.path.isdir(output + "/sample_counts_distrib"):
        #         os.mkdir(output + "/sample_counts_distrib")
        #     fig = ff.create_distplot([counts_persistent, counts_shell], ["persistent", "shell"],
        #                              colors=[COLORS_partitions["persistent"], COLORS_partitions["shell"]], bin_size=10)
        #     fig.update_layout(title=f"distribution of counts, sample {sample}")
        #     fig.add_vline(x=sample_th[sample], opacity=0.3, line_width=0.5)
        #     out_plotly.plot(fig, filename=output + "/sample_counts_distrib/" + sample + ".html", auto_open=False)
        # significant_families = {}
        # samples_dataset1 = set(samples_dataset1)
        # samples_dataset2 = set(samples_dataset2)
        # if filter_by_stats:
        #     for node in families:
        #         dataset1_fampresence = 0
        #         dataset2_fampresence = 0
        #         dataset1_famabsence = 0
        #         dataset2_famabsence = 0
        #         for sample in sample_index:
        #             if pangenome.samples[sample].gene_families_map_count[node.name] > sample_th[sample]:
        #                 if sample in samples_dataset1:
        #                     dataset1_fampresence += 1
        #                 elif sample in samples_dataset2:
        #                     dataset2_fampresence += 1
        #                 else:
        #                     import pdb;
        #                     pdb.set_trace()
        #             else:
        #                 if sample in samples_dataset1:
        #                     dataset1_famabsence += 1
        #                 elif sample in samples_dataset2:
        #                     dataset2_famabsence += 1
        #                 else:
        #                     import pdb;
        #                     pdb.set_trace()
        #         contingency_table = np.array(
        #             [[dataset1_fampresence, dataset2_fampresence], [dataset1_famabsence, dataset2_famabsence]])
        #         V = float("nan")
        #         try:
        #             oddsratio, pvalue = stats.fisher_exact(contingency_table)
        #             X2 = stats.chi2_contingency(contingency_table, correction=False)[0]
        #             n = np.sum(contingency_table)
        #             minDim = min(contingency_table.shape) - 1
        #             V = np.sqrt((X2 / n) / minDim)
        #         except:
        #             print(node.name + " " + str(contingency_table))
        #         significant_families[node] = [pvalue, oddsratio, V]
        # cor_p_values = \
        # multipletests(list([stats_info[0] for stats_info in significant_families.values()]), alpha=0.05, method='fdr_bh')[
        #     1].tolist()
        # for index_fam, fam in enumerate(significant_families.keys()):
        #     significant_families[fam] = significant_families[fam] + [cor_p_values[index_fam]]
        # families = [f for f, stats_info in significant_families.items()]
        persistent_fam = [fam for fam in families if fam.named_partition == "persistent"]
        shell_fam = [fam for fam in families if fam.named_partition == "shell"]
        cloud_fam = [fam for fam in families if fam.named_partition == "cloud"]
        # print([fam.name for fam in shell_fam+cloud_fam])
        # if filter_by_stats:
        data = []
        for sample in sample_index:
            data.append(
                [self.sample_counts_gene_families[fam][sample] if self.sample_counts_gene_families[fam][sample] > sample.th else 0 for fam in families])

        mat_mapping = np.matrix(data, dtype='bool')

        # mat_p_a = csc_matrix((data, (all_indexes,all_columns)), shape = (len(families),len(pangenome.samples)), dtype='float')
        # dist = pairwise_distances(mat_p_a)
        dist = pdist(mat_mapping, 'jaccard')
        hc = linkage(dist, 'ward')
        dendro = dendrogram(hc, no_plot=True)

        logging.getLogger().info("done with making the dendrogram to order the samples on the plot")
        order_samples = [index2sample[index] for index in dendro["leaves"]]

        #fig_samples = ff.create_dendrogram(mat_mapping, orientation='top', distfun=lambda x: pdist(x, metric='jaccard'),
        #                                   linkagefun=lambda x: hc, labels=[s.name for s in order_samples])
        #for i in range(len(fig_samples['data'])):
        #    fig_samples['data'][i]['yaxis'] = 'y2'

        sample_data = []
        text_data = []
        fam_order = []
        partitions_dict = defaultdict(list)
        shell_subs = set()
        ordered_nodes = []
        order_fam = []
        # for fam in families:
        #     sum_mapping = sum([pangenome.samples[sample].gene_families_map_count[fam.name] for sample in sample_index])
        #     partitions_dict[fam.partition].append([fam,sum_mapping])
        #     if fam.partition.startswith("S"):
        #         shell_subs.add(fam.partition)#number of elements will tell the number of subpartitions
        # ordored_nodes_p = sorted(partitions_dict["P"], key=lambda n:n[1], reverse=True)
        # ordored_nodes_c = sorted(partitions_dict["C"], key=lambda n:n[1], reverse=True)
        order_fam_persistent = []
        order_fam_shell = []
        order_fam_cloud = []

        if len(persistent_fam) > 0:
            data = []
            index2fam = {}
            for index, fam in enumerate(persistent_fam):
                index2fam[index] = fam
                data.append(
                    [1 if self.sample_counts_gene_families[fam][sample] > sample.th else 0 for sample in
                     order_samples])
            mat_mapping = np.matrix(data, dtype='bool')
            dist = pdist(mat_mapping, 'jaccard')
            hc = linkage(dist, 'ward')
            dendro = dendrogram(hc, no_plot=True)
            order_fam_persistent = [index2fam[index] for index in dendro["leaves"]]
            order_fam_persistent.reverse()
        if len(shell_fam) > 0:
            data = []
            index2fam = {}
            for index, fam in enumerate(shell_fam):
                index2fam[index] = fam
                data.append(
                    [1 if self.sample_counts_gene_families[fam][sample] > sample.th else 0 for sample
                     in order_samples])
            mat_mapping = np.matrix(data, dtype='bool')
            dist = pdist(mat_mapping, 'jaccard')
            hc = linkage(dist, 'ward')
            dendro = dendrogram(hc, no_plot=True)
            order_fam_shell = [index2fam[index] for index in dendro["leaves"]]
            order_fam_shell.reverse()
            order_fam_cloud = []
        if not nocloud:
            data = []
            index2fam = {}
            for index, fam in enumerate(cloud_fam):
                index2fam[index] = fam
                data.append(
                    [1 if self.sample_counts_gene_families[fam][sample] > sample.th else 0 for sample
                     in order_samples])
            mat_mapping  = np.matrix(data, dtype='bool')
            dist = pdist(mat_mapping, 'jaccard')
            hc = linkage(dist, 'ward')
            dendro = dendrogram(hc, no_plot=True)
            order_fam_cloud = [index2fam[index] for index in dendro["leaves"]]
            order_fam_cloud.reverse()

        sep_p = len(order_fam_persistent) - 0.5
        separators = [sep_p]
        shell_NA = None
        order_fam = order_fam_persistent + order_fam_shell + order_fam_cloud
        separators.append(separators[len(separators) - 1] + len(order_fam_shell))
        separators.append(separators[len(separators) - 1] + len(order_fam_cloud))

        logging.getLogger().info("Getting the gene name(s) and the number for each tile of the plot ...")
        to_print = []

        shapes = []
        for num, node in enumerate(order_fam):
            if self.results[node].pvalue < 0.05 :
                if self.results[node].oddsratio>1:
                    shapes.append(
                        dict(type='line', x0=-2, x1=-2, y0=num-0.5, y1=num+0.5, line=dict(dict(width=10, color="brown"))))
                    shapes.append(
                        dict(type='line', x0=len(order_samples)+2, x1=len(order_samples)+2, y0=num-0.5, y1=num+0.5, line=dict(dict(width=10, color="brown"))))
                    fam_order.append('<b><span style="color:brown">\u200c' + node.name + '</span>')
                elif self.results[node].oddsratio<1:
                    shapes.append(
                        dict(type='line', x0=-2, x1=-2, y0=num - 0.5, y1=num + 0.5,
                             line=dict(dict(width=10, color="limegreen"))))
                    shapes.append(
                        dict(type='line', x0=len(order_samples) + 2, x1=len(order_samples) + 2, y0=num - 0.5,
                             y1=num + 0.5, line=dict(dict(width=10, color="limegreen"))))
                    fam_order.append('<b><span style="color:limegreen">\u200c' + node.name + '</span>')
                else:
                    fam_order.append('<b>'+node.name+'</b>')
            else:
                fam_order.append('<b>'+node.name+'</b>')
            v = 0
            try:
                V = self.results[node].cramer_V
                # results[f].oddsratio
                # v = comparison.   significant_families[node][2] if significant_families[node][3] < 0.01 else 0
                # if V > 0.4:
                #    to_print.append(node.name)
            except:
                pass
            sample_data.append([self.sample_counts_gene_families[node][sample] if self.sample_counts_gene_families[node][sample] > sample.th else np.nan for sample in order_samples])
            text_data.append([f"{node.name}<br>\
                                {node.partition}<br>\
                                {np.format_float_scientific(self.sample_counts_gene_families[node][sample])}<br>\
                                pval={np.format_float_scientific(self.results[node].pvalue)}<br>\
                                OR={np.format_float_scientific(self.results[node].oddsratio)}<br>\
                                V={np.format_float_scientific(self.results[node].cramer_V)}<br>\
                                pval_c={np.format_float_scientific(self.results[node].corrected_pvalue)}" for sample in
                                order_samples])  # <br>{dataset[sample]}
        print("|".join(to_print))
        #import pdb;
        #pdb.set_trace()
        xaxis_values = ['<b><span style="color:brown">'+sample.name+'</span></b>' if sample.name in self.condition1_samples_or_orgs else '<b><span style="color:limegreen">'+sample.name+'</span></b>' for sample in order_samples]

        logging.getLogger().info("Done extracting names and numbers. Making the heatmap ...")

        heatmap = go.Heatmap(z=sample_data,
                             x=xaxis_values,
                             y=fam_order,
                             text=text_data,
                             zmin=0,
                             zmax=1,
                             autocolorscale=True)
        #heatmap['x'] = fig_samples['layout']['xaxis']['tickvals']
        #fig_samples.add_trace(heatmap)
        if len(shell_subs) > 1:
            if "S_" not in shell_subs:
                shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7], len(shell_subs))
            else:
                shell_color = cl.interp(cl.flipper()['seq']['9']['Greens'][1:7], len(shell_subs) - 1)

        sep_prec = 0
        for nb, sep in enumerate(separators):
            color = None
            if nb == 0:
                color = COLORS["persistent"]
            elif nb == (len(separators) - 1):
                color = COLORS["cloud"]
            elif len(shell_subs) > 1:
                if shell_NA is not None and nb == shell_NA:
                    color = COLORS["shell"]
                else:
                    color = shell_color.pop()
            else:
                color = COLORS["shell"]
            shapes.append(dict(type='line', x0=-1, x1=-1, y0=sep_prec, y1=sep, line=dict(dict(width=10, color=color))))
            shapes.append(dict(type='line', x0=len(order_samples)+1, x1=len(order_samples)+1, y0=sep_prec, y1=sep,
                               line=dict(dict(width=10, color=color))))
            shapes.append(
                dict(type='line', x0=-1, x1=len(order_samples)+1, y0=sep, y1=sep, line=dict(dict(width=1, color=color))))
            sep_prec = sep
        for index_order, sample in enumerate(order_samples):
            if sample.name in self.condition1_samples_or_orgs:
                shapes.append(dict(type='line', x0=index_order - 0.5, x1=index_order + 0.5, y0=-5, y1=-5, line=dict(dict(width=10, color="brown"))))
            elif sample.name in self.condition2_samples_or_orgs:
                shapes.append(dict(type='line', x0=index_order - 0.5, x1=index_order + 0.5, y0=-15, y1=-15, line=dict(dict(width=10, color="limegreen"))))

        layout = go.Layout(title="mapping matrix",
                           xaxis=go.layout.XAxis(ticktext=xaxis_values,
                                                 title='samples <b><span style="color:brown">'+self.condition1_name+'</span></b> vs <b><span style="color:limegreen">'+self.condition2_name+'</span></b>',
                                                 tickvals=xaxis_values,
                                                 automargin=True,
                                                 tickfont=dict(size=10)),
                           yaxis=go.layout.YAxis(ticktext=fam_order,
                                                 tickvals=fam_order,
                                                 title='gene families',
                                                 automargin=True,
                                                 tickfont=dict(size=10)),
                           shapes=shapes,
                           plot_bgcolor='#ffffff')
        logging.getLogger().info("Drawing the figure itself...")
        
        out_plotly.plot(go.Figure(data=[heatmap], layout=layout), filename=output + ".html", auto_open=False)
        logging.getLogger().info(f"Done with the sample tile plot : '{output + '.html'}' ")

def extract_modules_file(functional_modulesf_file):
    """ Extract the module associated to each families

    :param functional_modulesf_file: input path 
    :type functional_modulesf_file: str
    :return: a dictionary where families are keys and sets of modules are values 
    :rtype: dict
    """
    families_modules=defaultdict(set)
    with open(functional_modulesf_file,"r") as tsvfile:
        for i, line in enumerate(tsvfile):
            elements = [e.strip() for e in line.split("\t")]
            genomes_or_sample_name = elements[0]
            if i == 0:
                continue
            else:
                families_modules[elements[0]].add(elements[1])
    return(families_modules)

def extract_conditions_and_sample_read_files(pangenome : Pangenome, file : str):
    """ extract condition names and list of genomes for these two condition from matrice file 
    :param file: contains genome_names, condition1 and condition2 with 1 if the genome is present for condition1 or condition2 and 0 if the genome is absent
    :type str
    :return: a tuple where the 1st element is the list of sample or organisme for the condition 1, the 2nd is the same for the the condition 2, the 3rd is the name of the condition 1, the 4th is the name of the condition 2 and the last one is a dictionary where sample Id are keys and each value correspond the a list of the paths to reach the read files (index 0 for R1 and index 1 for R2)
    :rtype: tuple
    """
    #check_pangenome_info(pangenome, need_families=True, need_annotations=True)

    comparison = None
    with open(file,"r") as tsvfile:
        for i, line in enumerate(tsvfile):
            elements = [e.strip() for e in line.split('\t')]
            genomes_or_sample_name = elements[0]
            if i == 0:
                comparison = Comparison(pangenome = pangenome, condition1_name = elements[1], condition2_name = elements[2])
            else:
                try:
                    to_be_compared = pangenome.get_organism(genomes_or_sample_name)
                except KeyError as e:
                    to_be_compared = Sample(genomes_or_sample_name)
                    if len(elements) == 4:
                        to_be_compared.path_read1 = elements[3]
                    if len(elements) == 5:
                        to_be_compared.path_read1 = elements[3]
                        if elements[4] != "":
                            to_be_compared.path_read2 = elements[4]
                if (elements[1] == "1"):
                    comparison.add_sample_or_org(to_be_compared, condition = 1)
                if (elements[2] == "1"):
                    comparison.add_sample_or_org(to_be_compared, condition = 2)
                if (elements[1] == "1" and elements[2] == "1"):
                    raise logging.getLogger.warning(f"A Genome or a sample is present at same time for two conditions : '{genomes_or_sample_name}")

        logging.getLogger().debug("end of extraction step")
        return(comparison)

def launch(args):
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)

    check_pangenome_info(pangenome, need_families=True, need_annotations=True)
    if args.conditions_file is None or os.stat(args.conditions_file).st_size == 0:
      raise Exception(f"Your provided file is empty")

    comparison = extract_conditions_and_sample_read_files(pangenome, args.conditions_file)
    if args.import_count_matrix is not None:
        comparison.import_count_matrix(count_matrix_file=args.import_count_matrix)
    if len(comparison.get_all_samples_having_counts(reversed=True)) > 0:
        map_on_graphs(pangenome, comparison, output_dir = args.output, min_ident=args.min_identity, cpu= args.cpu)
        with open(args.output+"/count_matrix_"+comparison.name+".tsv","w") as count_matrix:
            samples = list([s for s in comparison.get_all_samples()])
            count_matrix.write("\t".join(["gene_families"]+[s.name for s in samples])+"\n")
            for gf in pangenome.gene_families:

                try:
                    count_matrix.write("\t".join([str(gf.name)]+[str(comparison.sample_counts_gene_families[gf][s]) for s in samples])+"\n")
                except:
                    import pdb;
                    pdb.set_trace()
    #TODO retreive this information from .h5 file
    functional_modules=None
    if args.import_functional_modules is not None:
        functional_modules = extract_modules_file(args.import_functional_modules)
    logging.getLogger().debug("start of perform_compararison step")
    comparison.perform_comparisons(functional_modules = functional_modules,
                                   min_cov_persistent = args.min_cov_persistent,
                                   th_ratio_persistent_mean_coverage = args.th_ratio_persistent_mean_coverage,
                                   dir = args.output)

    comparison.draw_tileplot_samples(args.output+"/tile_plot_samples_"+comparison.name, nocloud=args.plot_nocloud, filter_by_stats=True)
    logging.getLogger().debug("end of performCompararison step")

def subparser(sub_parser):
    """
    Parser arguments specific to compare command

    :param sub_parser : sub_parser for align command
    :type sub_parser : argparse._SubParsersAction

    :return : parser arguments for align command
    :rtype : argparse.ArgumentParser
    """
    parser = sub_parser.add_parser("compare", formatter_class=argparse.RawTextHelpFormatter)
    parser_compare(parser)
    return parser

def parser_compare(parser):

    required = parser.add_argument_group(title = "Required arguments",
                                         description = "All of the following arguments are required :")
    required.add_argument('-p', '--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str, help="Output directory where the file(s) will be written")
    required.add_argument('-cf', '--conditions_file', required=True, type=str, help='''The elements included in the two conditions (either genomes or samples, mixed or not) and if necessary R1 and R2 reads paths of the samples to map them against the pangenome. Here is a toy example file (header required and elements must be tab-separated):
    sample_or_organism  condition1 condition2 R1                 R2
    sample1             1          0          sample1_1.fastq.gz sample1_2.fastq.gz
    sample2             0          1          sample2_1.fastq.gz sample2_2.fastq.gz
    sample3             1          0          sample3_1.fastq.gz sample3_2.fastq.gz
    organism_1          1          0
    organism_2          0          1''')

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--import_count_matrix', default=None, required = False, type = str, help ='''Allow to import a previously computed mapping of samples against the pangenome families. Here is a toy example file (header required and elements must be tab-separated):
    families  sample1 sample2 sample3
    families1 1.2     1.5     0.6
    families2 4.1     2.1     0.8
    families3 2.8     0       0.4''')
    optional.add_argument('--import_functional_modules', default=None, required = False, type = str, help = "Allow to import the list of module associated with each families. This correspond to the \"functional_modules.tsv\" file generate by the sub command \"module\"")
    optional.add_argument('--min_cov_persistent',required=False, type=restricted_float, default=0.85, help = "Specify the ratio of persistent genome to be covered to include a sample in the comparison")
    optional.add_argument('--th_ratio_persistent_mean_coverage', required=False, type=restricted_float, default=0.05, help = "The percentile of the persistent genome distribution which coverage above means presence and below means absence")
    optional.add_argument("--plot_nocloud", required=False, default=False, action="store_true",
                          help="Do not draw the cloud in the tile plot")
    optional.add_argument("--min_identity", required=False, type=restricted_float, default=0.95, help = "Accept each read only if the alignment identity is above this minimum")
    return parser

if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_compare(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())