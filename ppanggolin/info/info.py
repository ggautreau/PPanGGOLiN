#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
from pathlib import Path
import logging

# installed libraries
import tables
import yaml

# local libraries
from ppanggolin.formats import read_info, read_parameters,  fix_partitioned


def read_status(h5f: tables.File):
    status_group = h5f.root.status

    status_to_print = {
        "Genomes_Annotated": status_group._v_attrs.genomesAnnotated,
        "Genes_Clustered": status_group._v_attrs.genesClustered,
        "Genes_with_Sequences": status_group._v_attrs.geneSequences,
        "Gene_Families_with_Sequences": status_group._v_attrs.geneFamilySequences,
        "Neighbors_Graph": status_group._v_attrs.NeighborsGraph,
        "Pangenome_Partitioned": status_group._v_attrs.Partitioned,
        "RGP_Predicted": status_group._v_attrs.predictedRGP,
        "Spots_Predicted": status_group._v_attrs.predictedRGP, # Please confirm if this should be different from "RGP Predicted"
        "Modules_Predicted": status_group._v_attrs.modules
    }
    status_to_print = {key:bool(val) for key, val in status_to_print.items()}

    status_to_print["PPanGGOLiN_Version"] = str(status_group._v_attrs.version)

    yaml_output = yaml.dump({"Status":status_to_print}, default_flow_style=False, sort_keys=False, indent=4)
    print(yaml_output)

def read_metadata(h5f):
    status_group = h5f.root.status
    if hasattr(status_group._v_attrs, "metadata") and status_group._v_attrs.metadata:
        metastatus = status_group.metastatus
        metasources = status_group.metasources
        print("Metadata: ")
        for attr in metastatus._v_attrs._f_list():
            print(f"    - {attr} : {', '.join(metasources._v_attrs[attr])}")
    else:
        logging.getLogger("PPanGGOLiN").warning("There is not any metadata in the pangenome")


def print_info(pangenome: str, status: bool = False, content: bool = False, parameters: bool = False,
               metadata: bool = False):
    """
    Main function to return information about pangenome

    :param pangenome: Pangenome file
    :param status: Get pangenome status
    :param content: Get pangenome content
    :param parameters: Get pangenome parameters
    """
    fix_partitioned(pangenome)
    if status or content or parameters or metadata:
        h5f = tables.open_file(pangenome, "r+")
        if status:
            read_status(h5f)
        if content:
            read_info(h5f)
        if parameters:
            read_parameters(h5f)
        if metadata:
            read_metadata(h5f)
        h5f.close()
    else:
        raise ValueError("Please select what information you want by using --parameters, --content or --status")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    print_info(args.pangenome, args.status, args.content, args.parameters, args.metadata)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for info command

    :return : parser arguments for info command
    """
    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_info(parser)
    return parser


def parser_info(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of info command

    :param parser: parser for info argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="The following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=Path, help="The pangenome .h5 file")

    options = parser.add_argument_group(title="optional arguments")
    options.add_argument("--parameters", required=False, action="store_true",
                         help="Shows the parameters used (or computed) for each step of the pangenome generation")
    options.add_argument("--content", required=False, action="store_true",
                         help="Shows detailled informations about the pangenome's content")
    options.add_argument("--status", required=False, action="store_true",
                         help="Shows informations about the statuses of the different elements of the pangenome "
                              "(what has been computed, or not)")
    options.add_argument("--metadata", required=False, action="store_true",
                         help="Shows which metadata are saved in the pangenome")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_info(main_parser)
    launch(main_parser.parse_args())
