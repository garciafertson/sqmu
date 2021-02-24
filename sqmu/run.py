#!/usr/bin/env python
import subprocess
import os
import sys
from sqmu.modules import run_kmod_taxlev, run_kmod_taxlist, prep_kmfiles
from sqmu.modules_keggmapper import run_kmod_kmap
from sqmu.network import build_network
from sqmu.graftm import graftm_abundance


class Run:
    def __init__(self, args):
        self.args = args

    # workflow to build table
    def build_tax_module(self):
        # save parameter in variables
        sqzm_folder = self.args.sqzm_folder
        modules = self.args.modules
        tax_file = self.args.taxonomy_list_file
        bin_list_file = self.args.bin_list_file
        out_table = self.args.out_table
        custom_mod = self.args.custom_mod
        if not self.args.prepare_files_keggmap:
            if len(self.args.keggmapper) == 0:
                if self.args.taxlev:
                    run_kmod_taxlev(sqzm_folder, modules,custom_mod, self.args.taxlev,out_table )
                else:
                    run_kmod_taxlist(sqzm_folder, modules, custom_mod, tax_file,
                                     bin_list_file, out_table)
            else:
                km_files = self.args.keggmapper
                run_kmod_kmap(km_files,modules,out_table)
        else:
            prep_kmfiles(sqzm_folder, tax_file, bin_list_file)

    # workflow to buil network files

    def build_ko_network(self):
        sqzm_folder = self.args.sqzm_folder
        kolist_file = self.args.kolist_file
        taxlevel = self.args.taxlevel
        out=self.args.output
        build_network(sqzm_folder, kolist_file, taxlevel, out)

    # workflow to add abundance data into graftM results
    def graftm_ab(self):
        sqzm_folder = self.args.sqzm_folder
        graftm_result = self.args.graftm_result
        abundance = self.args.abundance
        out = self.args.out
        graftm_abundance(sqzm_folder, graftm_result, abundance, out)

    def main(self):
        if self.args.subparser_name == 'tax_module':
            self.build_tax_module()
            print('MODULE_TABLE done')
        elif self.args.subparser_name == 'tax_network':
            self.build_ko_network()
            print('KO_NETWORK done')
        elif self.args.subparser_name == 'graftm_abundance':
            self.graftm_ab()
            print('GRAFTM SQM ABUNDANCE done')
