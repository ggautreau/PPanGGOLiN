#!/usr/bin/env python3
#coding: utf8

#default libraries
import logging
from collections.abc import Iterable
from scipy.stats import ttest_rel, ttest_ind
from statsmodels.stats.multitest import multipletests

#local libraries
from ppanggolin.sample import Sample

class Dataset:
    """This represents a 
    """

    def __init__(self, dataset, samples_dataset = set()):
        """Constructor method
        :param ID: The sample ID
        :type name: str
        """
        self.ID = dataset
        self.samples_dataset = samples_dataset

    @property
    def gene_families_map_count_mean(self):
        try:
            return self._gene_families_map_count_meanGetter
        except AttributeError:#in that case the attribute has not been computed
            self._mkgene_families_map_count_meanGetter()
            return self._gene_families_map_count_meanGetter#return what was expected

    def _mkgene_families_map_count_meanGetter(self):
        self._gene_families_map_count_meanGetter = {}
        for fam, map_count_mean in self._yield_gene_families_map_count_mean():
            self._gene_families_map_count_meanGetter[fam] = map_count_mean
        
    def _yield_gene_families_map_count_mean(self):
        fams = list(self.samples_dataset)[0].gene_families_map_count.keys()
        for fam in fams:
            sum_count = 0
            for sample in self.samples_dataset:
                sum_count += sample.gene_families_map_count[fam]
            yield(tuple([fam, sum_count / len(self.samples_dataset)]))

class Comparison:
    """This represents a 
    """

    def __init__(self, IDcomparison, dataset1 = [], dataset2 = [], paired = False, corrected_p_values = True):
        """Constructor method
        :param IDcomparison: The comparison ID
        :type name: str
        """
        self.ID = IDcomparison
        self.dataset1 = dataset1
        self.dataset2 = dataset2
        self.paired = paired
        self.corrected_p_values = corrected_p_values

    @property
    def gene_families_map_count_mean_FC(self):
        try:
            return self._gene_families_map_count_mean_FCGetter
        except AttributeError:#in that case the attribute has not been computed
            self._mk_gene_families_map_count_mean_FCGetter()#return what was expected
            return self._gene_families_map_count_mean_FCGetter

    def _mk_gene_families_map_count_mean_FCGetter(self):
        self._gene_families_map_count_mean_FCGetter = {}
        for fam, FC in self._yield_gene_families_map_count_mean_FC():
            self._gene_families_map_count_mean_FCGetter[fam] = FC

    def _yield_gene_families_map_count_mean_FC(self):
        fams = list(self.dataset1.samples_dataset)[0].gene_families_map_count.keys()
        for fam in fams:
            FC = 0
            try :
                FC = self.dataset1.gene_families_map_count_mean[fam] / self.dataset2.gene_families_map_count_mean[fam] 
            except ZeroDivisionError:
                pass
            yield(tuple([fam, FC]))

    @property
    def gene_families_map_count_p_values(self):
        try:
            return self._gene_families_map_count_p_valuesGetter
        except AttributeError:#in that case the attribute has not been computed
            self._mk_gene_families_map_count_p_valuesGetter()#return what was expected
            return self._gene_families_map_count_p_valuesGetter

    def _mk_gene_families_map_count_p_valuesGetter(self):
        self._gene_families_map_count_p_valuesGetter = {}
        fams = list(self.dataset1.samples_dataset)[0].gene_families_map_count.keys()
        for fam in fams:
            if self.paired:
                self._gene_families_map_count_p_valuesGetter[fam] = ttest_rel([sample.gene_families_map_count[fam] for sample in self.dataset1.samples_dataset],
                                                                              [sample.gene_families_map_count[fam] for sample in self.dataset2.samples_dataset]).pvalue
            else:
                self._gene_families_map_count_p_valuesGetter[fam] = ttest_ind([sample.gene_families_map_count[fam] for sample in self.dataset1.samples_dataset],
                                                                              [sample.gene_families_map_count[fam] for sample in self.dataset2.samples_dataset]).pvalue
            
            print([sample.gene_families_map_count[fam] for sample in self.dataset1.samples_dataset])
            exit
        if self.corrected_p_values:
            cor_p_values = multipletests(self._gene_families_map_count_p_valuesGetter.values(), alpha=0.05, method='fdr_bh')
            self._gene_families_map_count_p_valuesGetter = dict(zip(self._gene_families_map_count_p_valuesGetter.keys(), cor_p_values))