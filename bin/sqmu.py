#!/usr/bin/env python
################################################
# Squeeze meta some exploratory utilities
################################################
import logging
import os
import sys
import argparse

__author__ = "Fernando Garcia Guevara"
__credits__ = "Fernando Garcia Guevara"
__email__ = "garciafertson near gmail.com"
__status__ = "Development"

# this line get the path where the program is installed
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')]+sys.path
# import module for running program
import sqmu
from sqmu.run import Run

def print_header():
    print """ Squeeze Meta Utilities """


def phelp():
    print"""
        SQMU
        This program performs a simple exploratoty analysis on Squeezemeta v1.3
        results. For now there are only two types of analysis and all
        of them are focused on KO and taxonomy SQM results folder
        resylts.
        For more information please type:

        sqmu.py tax_module -h

        sqmu.py tax_network -h

        sqmu.py graftm_abundance -h
    """


# the next section add parsers and subparser to the program
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='sqmu v%s' % sqmu.__version__)
    subparser = parser.add_subparsers(help="sub-comand help:", dest='subparser_name')

    # PARSER 1 taxa/bin modules parser
    taxmod_parser = subparser.add_parser('tax_module',
                                         description='Build tsv table of modules against taxonomy',
                                         help='Input squeezemeta results table and specify options',
                                         epilog='''This option creates a tsv with presecene/absence values for the specifed kegg modules.
    You need internete connection ass this tool uses KEGGREST api to download module information ''')

    input_options = taxmod_parser.add_argument_group('input options')
    input_options.add_argument("-f", '--sqzm_folder',
                               metavar='FOLDER PATH',
                               help="Path to Squeezemeta output folder",
                               required=True)
    input_options.add_argument("-m", '--modules',
                               metavar="FILE",
                               help="File with list of Kegg modules to be included in the table",
                               required=True)
    input_options.add_argument('--taxonomy_list_file',
                               metavar='FILE',
                               help="(optional) File with list of NCBI taxonomy levels to be included in table",
                               default="")
    input_options.add_argument('--bin_list_file',
                               metavar='FILE',
                               help="(optional) File with list of bin names to be included in table",
                               default="")
    input_options.add_argument('--custom_mod',
                               metavar="FILE",
                               help="(optional) File with custome modules, the format consist of a plain text,\
        a module starts with > and then the descriptor or module name, then each line \
        specify the module steps, conformed by KO separated by comma, until the next > ",
                               default="")
    input_options.add_argument('--keggmapper',
                               nargs="*",
                               metavar='FILE LIST',
                               help="(optional)This option overrides previous options. Space separated list of files with the results from the KeggMapper WEB tool for each organism/bin to be included in the final table, as there is no automatic way (for now) to make queries you should do it         manually and save the results in a plain text file",
                               default=[])
    input_options.add_argument("--taxlev",
                               metavar='STR',
                               help="deepest taxonomic level to construct module table\
                            p_;c_;o_;f_;g_;s_ this option overrides \
                            --taxonomy_file_list and creates a file with taxonomy list")

    output_options = taxmod_parser.add_argument_group("Output options")
    output_options.add_argument("-o", "--out_table",
                                metavar='STR',
                                default="output_module_table.tsv",
                                help="set output file name prefix")

    pf_options = taxmod_parser.add_argument_group("prepare files for Keggmaper")
    pf_options.add_argument("--prepare_files_keggmap",
                            action='store_true',
                            help="if specified it produces a list of files ready to be uploaded to the \
        KEGG MAPPER tool. You should load each resulting file into the tool \
        (https://www.genome.jp/kegg/tool/map_pathway.html). In the resulting page select \
        the Module tab and view inlcuding any incomplete, and save results in a plain text file")

    # PARSER 2 Build taxonomy and KO network files
    network_subparser = subparser.add_parser('tax_network',
                                             description="Build the Network files in paired format ready to be imported into Cytoscape program",
                                             help="This option recieves a list of KO and a taxonomic level",
                                             epilog=''' This option recieves a list ofKO and a taxonomic level, then searches in the orftable\
        SQM results folder and builds a network where KO and taxonomic categories are represented as nodes\
        and  edges represent the abundance found for a KO and taxonomic categorie pair''')

    netinput_options = network_subparser.add_argument_group("input files for network")
    netinput_options.add_argument("-f", "--sqzm_folder",
                                  metavar="FOLDER PATH",
                                  help="Path to Squeezemeta output folder",
                                  required=True)
    netinput_options.add_argument("-k", "--kolist_file",
                                  metavar="FILE",
                                  help="File with list of KO numbers to be included in the network",
                                  required=True)
    netinput_options.add_argument("-l", "--taxlevel",
                                  metavar="STR",
                                  help="Deepest taxonomic level to include in the network p_;c_;o_;f_;g_;s_ (default c_)",
                                  default="c_")
    netoutput_options=network_subparser.add_argument_group("output options for network")
    netoutput_options.add_argument("-o","--output",
                                    metavar="STR",
                                    help="Prefix for output files (default network)",
                                    default="network")

    # PARSER3 add tpm abundance to graftM results
    graftm_subparser = subparser.add_parser('graftm_abundance',
                                            description="add tpm abundance values to graftM results",
                                            help="This option runs grafm on sqzm predicted genes and add coverage values",
                                            epilog='''Recieves the graftm output ran on the results file
        03.$sqm_project_name.fna, the szqm output folder and add abundance values to the
        graftm output. For the moment only works with SQZM runs with onlu one sample\
        ''')

    graftm_options = graftm_subparser.add_argument_group("input graftm abundace")
    graftm_options.add_argument("-f", "--sqzm_folder",
                                metavar="FOLDER PATH",
                                help='Path to Squeezemeta output folder',
                                required=True)
    graftm_options.add_argument("--graftm_result",
                                metavar="FOLDER PATH",
                                help='Path to GraftM output folder',
                                required=True)
    graftm_output = graftm_subparser.add_argument_group("output options")
    graftm_output.add_argument("--abundance",
                               metavar="STR",
                               help="select abundance value, possible values are: Reads, RPKM and TPM (default TPM)",
                               default="TPM")
    graftm_output.add_argument("--out",
                               metavar="STR",
                               help="output filename prefix",
                               default="graftm_abundance_table")

    # check whether --help is needed
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    # call Run module passing the arguments in here
    else:
        args = parser.parse_args()
        Run(args).main()
