#!/usr/bin/env python3
#coding: utf8

#default libraries
import logging
from collections.abc import Iterable

#local libraries
#from ppanggolin.geneFamily import GeneFamily

class Sample:
    """This represents a sample (reads) mapped over the pangenome
    """

    def __init__(self, ID):
        """Constructor method
        :param ID: The sample ID
        :type name: str
        """
        self.ID = ID
        self.gene_families_map_count = dict()