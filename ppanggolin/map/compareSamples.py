#!/usr/bin/env python3
#coding:utf-8

#default libraries
import logging
import tempfile
import subprocess
import argparse
from collections import defaultdict
import pysam
import sys

import scipy
import statsmodels
import plotly.graph_objs as go
import plotly.express as px
import plotly.offline as out_plotly
import math

#local libraries
from ppanggolin.formats import checkPangenomeInfo, writePangenome
from ppanggolin.utils import mkOutdir, read_compressed_or_not
from ppanggolin.pangenome import Pangenome
from ppanggolin.sample import Sample
from ppanggolin.compare import Dataset, Comparison

def compare(pangenome, dataset1, dataset2, samples_dataset1, samples_dataset2, paired = False, correct_p_values = True):
    # if pangenome.status["geneFamilySequences"] not in ["inFile","Loaded","Computed"]:
    #     print(pangenome.status["geneFamilySequences"])
    #     raise Exception("Cannot use this function as your pangenome does not have gene families representatives associated to it. For now this works only if the clustering is realised by PPanGGOLiN.")
    checkPangenomeInfo(pangenome, needFamilies=True, needSamples=True, needComparisons=True)
    # print(pangenome.samples)
    # print(samples_dataset1)
    # print(dataset1)
    # print([pangenome.samples[s] for s in samples_dataset1])
    # print(samples_dataset2)
    # print(dataset2)
    # print([pangenome.samples[s] for s in samples_dataset2])
    print(pangenome.samples)
    ds1 = Dataset(dataset1)
    ds1.samples_dataset = set([pangenome.samples[s] for s in samples_dataset1])
    ds2 = Dataset(dataset2)
    ds2.samples_dataset = set([pangenome.samples[s] for s in samples_dataset2])
    pangenome.addDataset(ds1)
    pangenome.addDataset(ds2)
    print(correct_p_values)
    c = Comparison(dataset1 + "_versus_" + dataset2, ds1, ds2)
    c.paired = paired
    c.corrected_p_values = correct_p_values
    pangenome.addComparison(c)

    pangenome.status["comparisons"] = "Computed"
    return(c)

def draw_volcano_plot(pangenome, comparison, output, FC_th = 2, p_values_th = 0.05):

    COLORS_partitions = {"shell": "#00D860", "persistent":"#F7A507", "cloud":"#79DEFF", "undefined":"#828282"}
    data_plot = []
    #if pangenome:#if partitioned
    corrected = 'C' if comparison.corrected_p_values else 'Not c'
    FCs = defaultdict(list)
    p_values = defaultdict(list)
    significants = defaultdict(list)
    
    texts = defaultdict(list)
    fam_list = []
    for fam in pangenome.geneFamilies:
        fc = comparison.gene_families_map_count_mean_FC[fam.name]
        pv = comparison.gene_families_map_count_p_values[fam.name]
        FCs[fam.namedPartition].append(fc)
        p_values[fam.namedPartition].append(pv)
        significants[fam.namedPartition].append("red" if ((fc >= FC_th or fc <= 1/FC_th) and pv <= p_values_th) else "DarkSlateGrey")
        texts[fam.namedPartition].append(f"""<b>{fam.name} </b><br> Number of genes : {len(fam.genes)}<br>Fold-change : {fc} <br>{corrected}corrected p_values : {pv}""")
        fam_list.append(f"""{fam.name}\t{fam.namedPartition}\t{len(fam.genes)}\t{fc}\t{pv}""")
        #print(FCs)
    for partition in ["undefined", "cloud", "shell", "persistent"]:
        #print(FCs[partition])
        #print(p_values[partition])
        data_plot.append(go.Scatter(x = FCs[partition],
                                    y = p_values[partition],
                                    mode = 'markers',
                                    marker = dict(opacity=0.9,line=dict(width=0.5, color=significants[partition])),
                                    name = partition,
                                    marker_color = COLORS_partitions[partition],
                                    text = texts[partition]))
    layout =  go.Layout(title = f"Volcano plot of the '{comparison.ID}' comparison",
                        xaxis = dict(title=f"Fold change (log2 scale)<br>{comparison.dataset2.ID} <= ENRICHED => {comparison.dataset1.ID}", type="log", dtick = math.log(2)),
                        yaxis = dict(title=corrected+"orrected p_values (reversed log10 scale)", type="log", dtick = math.log(10),autorange = "reversed"),
                        plot_bgcolor = '#ffffff')
    fig = go.Figure(data = data_plot, layout = layout)
    fig.add_hline(y=p_values_th, line_dash="dot")
    fig.add_vline(x=FC_th, line_dash="dot")
    fig.add_vline(x=1/FC_th, line_dash="dot")

    # shapes = [dict(type='line', x0=math.log2(FC_th), x1=math.log2(FC_th), y0=0, y1=max(list(comparison.gene_families_map_count_p_values.values())), line = dict(dict(width=1, dash='dashdot', color="grey"))),
    #                               dict(type='line', x0=math.log2(1/FC_th), x1=math.log2(1/FC_th), y0=0, y1=math.log10(max(list(comparison.gene_families_map_count_p_values.values()))), line = dict(dict(width=1, dash='dashdot', color="grey"))),
    #                               dict(type='line', x0=min(list(comparison.gene_families_map_count_mean_FC.values())), x1=max(list(comparison.gene_families_map_count_mean_FC.values())), y0=p_values_th, y1=p_values_th, line = dict(dict(width=5, dash='dashdot', color="grey")))],
                        
    out_plotly.plot(fig, filename = output+"/volcanoplot_"+comparison.ID+".html", auto_open=False)
    logging.getLogger().info(f"Done drawing the Volcano plot : '{output}/volcanoplot_{comparison.ID}.html'")
    with open(output+"/output_comparison.txt",'w') as output_comparison_file:
        output_comparison_file.write("\n".join(fam_list))

def launch(args):
    mkOutdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    if args.dataset1 is not None and args.dataset2 is not None and args.samples_dataset1 is not None and args.samples_dataset1 is not None:
        print(args.no_pvalues_correction)
        print(args.paired)
        c = compare(pangenome = pangenome, dataset1 = args.dataset1, dataset2 = args.dataset2, samples_dataset1 = args.samples_dataset1, samples_dataset2 = args.samples_dataset2, paired = args.paired, correct_p_values = args.no_pvalues_correction)
        writePangenome(pangenome, pangenome.file, args.force, show_bar = args.show_prog_bars)
        draw_volcano_plot(pangenome, c, args.output)

def compareSubparser(subparser):
    parser = subparser.add_parser("compare", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title = "Required arguments", description = "All of the following arguments are required :")
    onereq = parser.add_argument_group(title = "Input file", description = "One of the following argument is required :")
    onereq.add_argument('--dataset1', required = True, type = str, help = "name of the first dataset")
    onereq.add_argument('--dataset2', required = True, type = str, help = "name of the second dataset")
    onereq.add_argument('--samples_dataset1', required = True, nargs = "+", type = str, help = "list samples of the first dataset")
    onereq.add_argument('--samples_dataset2', required = True, nargs = "+", type = str, help = "list samples of the second dataset")
    onereq.add_argument('--paired', required = False, type = str, help = "paired datasets")
    onereq.add_argument('--no_pvalues_correction', default = True, action='store_false', help = "no p_values correction")

    required.add_argument('-p','--pangenome',  required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o','--output', required=True, type=str, help="Output directory where the file(s) will be written")
    return parser