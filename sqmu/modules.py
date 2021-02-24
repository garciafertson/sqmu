# Create an structure for module, with a list of steps and the
import os
import pickle
import re
import time
import numpy as np
from keggrest import keggrest as kegg
from sqmu.sqmfiles import orftable, bincontig


def module_keggrest(mod):
    try:
        entry = kegg.KEGGget(mod)
        time.sleep(2.5)
    except:
        print("module was not downloaded: " + mod)
    try:
        os.mkdir("tmp")
        print("tmp folder created")
    except:
        print("saving module in tmp folder")
    tmpdir = os.path.join(os.path.realpath("."), "tmp")
    outfile = os.path.join(tmpdir, mod+".pkl")
    with open(outfile, "w") as f:
        pickle.dump(entry, f)
    return(entry)


class kmodule:
    def __init__(self, k_id, name):
        self.name = name
        self.id = k_id
        self.kolist = []
        self.steps = []

    def read(self, path):
        reko = re.compile(r'K\d{5}')
        try:
            entry = os.join("tmp", self.id+".pkl")
            with open(entry, "r") as f:
                mod = pickle.load(f)
        except:
            mod = module_keggrest(self.id)
        definition = mod["DEFINITION"]
        orth = mod["ORTHOLOGY"]
        self.kolist = reko.findall(definition[0])
        steps = definition[0].split(" ")
        for step in steps:
            step_ko = set(reko.findall(step))
            self.steps.append(step_ko)

    def read_custome(self, line):
        line = line.rstrip()
        kos = line.split(",")
        self.kolist.extend(kos)
        self.steps.append(set(kos))


class tax:
    def __init__(self, name):
        self.name = name
        self.ko = set()

    def add_ko(self,ko):
        self.ko.add(ko)

    def get_ko_list(self, orf_table):
        next(orf_table)
        for line in orf_table:
            if self.name in line:
                ko = line.split("\t")[9]  # kegg id
                ko = ko.replace("*", "")
                if len(ko) > 0:
                    self.ko.add(ko)

    def print_ko(self):
        outname = self.name+".kolist"
        i = 0
        with open(outname, "w") as f:
            for ko in self.ko:
                f.write("%s\t%s\n" %(i, ko))
                i += 1


class mag:
    def __init__(self, name):
        self.name = name
        self.contigs = set()
        self.ko = set()

    def get_ko_list(self, bin_contig, orf_table):
        next(bin_contig)
        next(orf_table)
        for line in bin_contig:
            contig = line.split("\t")[0]
            biname = line.split("\t")[6]
            #print(biname, self.name)
            if self.name in biname:
                self.contigs.add(contig)
        # print(self.contigs)
        for line in orf_table:
            ko = line.split("\t")[9]
            ko = ko.replace("*", "")
            contig = line.split("\t")[1]
            if contig in self.contigs:
                self.ko.add(ko)

    def print_ko(self):
        outname = self.name+".kolist"
        i = 0
        with open(outname, "w") as f:
            for ko in self.ko:
                f.write("%s\t%s\n" %(i, ko))
                i += 1


def build_table(bin_dict, mod_dict, filename):
    with open(filename, "w") as o:
        header = ",".join(mod_dict.keys())
        o.write("bin_id,%s\n" % header)
        for bin_id in bin_dict.keys():
            print(bin_id, "\n")
            presence = np.zeros(len(mod_dict.keys()))
            bin_ko = bin_dict[bin_id].ko
            for i, module in enumerate(mod_dict.keys()):
                total_steps = float(len(mod_dict[module].steps))
                steps_found = float(0)
                # print(module)  # comment
                koinbin = bin_ko.intersection(mod_dict[module].kolist)
                koinbin = "+".join(koinbin)
                print(koinbin)
                print(mod_dict[module].steps)
                for step in mod_dict[module].steps:
                    intersect = bin_ko.intersection(step)
                    if len(intersect) > 0:
                        steps_found += 1  # step is considered to be found if at least 1 ko
                # categories
                fracc = steps_found/total_steps
                if fracc < 0.33:  # absent 0 < 0.5 fraction of steps found
                    presence[i] = 0
                elif fracc < 0.66:  # present 1 >= 0.5 and < than 0.7
                    presence[i] = 1
                elif fracc < 1:  # present 2 >= 0.7 and < than 1
                    presence[i] = 2
                else:  # complete 3 >= 1 fraction
                    presence[i] = 3
            values = ",".join([str(x) for x in presence])
            o.write("%s,%s\n" % (bin_id, values))


def run_kmod_taxlist(sqzm, mod, cust, taxlist, bin, out):
    print("searching presence of kegg modules")
    tax_dict = {}
    bin_dict = {}
    modules_dict = {}
    mod_list = open(mod, "r")
    # fill in dictionary with modules in list
    for line in mod_list:
        line = line.rstrip()
        mod_id = line.split("\t")[0]
        mod_name = line.split("\t")[1]
        modules_dict[mod_name] = kmodule(mod_id, mod_name)
        modules_dict[mod_name].read(mod_id)
    if len(cust) > 0:
        with open(cust) as custmod:
            for line in custmod:
                line = line.rstrip()
                if re.search(r"^>", line):
                    line.replace(">", "")
                    mod_id = line.split(",")[0]
                    mod_name = line.split(",")[1]
                    modules_dict[mod_name] = kmodule(mod_id, mod_name)
                else:
                    modules_dict[mod_name].read_custome(line)
    # read tax list
    if len(taxlist) > 0:
        tax_list = open(taxlist, "r")
        orftbl = orftable(sqzm)
        for line in tax_list:
            orf = open(orftbl, "r")
            line = line.strip()
            tax_dict[line] = tax(line)
            tax_dict[line].get_ko_list(orf)
            orf.close()
        tax_list.close()
        out_tax = out+"_tax.tsv"
        build_table(tax_dict, modules_dict, out_tax)

    if len(bin) > 0:
        bin_list = open(bin, "r")
        orftbl = orftable(sqzm)
        bintbl = bincontig(sqzm)
        for line in bin_list:
            orf = open(orftbl, "r")
            bin = open(bintbl, "r")
            line = line.strip()
            bin_dict[line] = mag(line)
            bin_dict[line].get_ko_list(bin, orf)
            orf.close()
            bin.close()
        bin_list.close()
        out_bin = out+"_bin.tsv"
        build_table(bin_dict, modules_dict, out_bin)


def run_kmod_taxlev(sqzm, mod, cust, taxlev, out):
    print("searching keggmodules in specified tax level")
    tax_dict = {}
    modules_dict = {}
    mod_list = open(mod, "r")
    # fill in dictionary with modules in list
    for line in mod_list:
        line = line.rstrip()
        mod_id = line.split("\t")[0]
        mod_name = line.split("\t")[1]
        modules_dict[mod_name] = kmodule(mod_id, mod_name)
        modules_dict[mod_name].read(mod_id)
    if len(cust) > 0:
        with open(cust) as custmod:
            for line in custmod:
                line = line.rstrip()
                if re.search(r"^>", line):
                    line.replace(">", "")
                    mod_id = line.split(",")[0]
                    mod_name = line.split(",")[1]
                    modules_dict[mod_name] = kmodule(mod_id, mod_name)
                else:

                    modules_dict[mod_name].read_custome(line)
    # read tax list
    orftbl = orftable(sqzm)
    orf = open(orftbl, "r")
    next(orf)
    next(orf)
    allkos = set()
    for key in modules_dict.keys():
        for ko in modules_dict[key].kolist:
            allkos.add(ko)
    print(allkos)
    for line in orf:
        line = line.strip()
        orf_annot = line.split("\t")
        taxonomy = orf_annot[8].split(";")
        ko = orf_annot[9]
        ko = ko.replace("*", "")
        if ko in allkos:
            if taxlev not in orf_annot[8]:
                if len(orf_annot[8]) == 0:
                    org = "Unclassified"
                else:
                    org = taxonomy[-1]  # select last taxonomy level
            else:
                for level in taxonomy:  # select specified taxonomy level
                    if taxlev in level:
                        org = level
            print(org,taxlev, tax_dict.keys())
            if org not in tax_dict.keys():
                tax_dict[org]=tax(org)
            tax_dict[org].ko.add(ko)

    out_tax = out+"_tax.tsv"
    build_table(tax_dict, modules_dict, out_tax)
    with open("taxlist.txt", "w") as t:
        for taxa in tax_dict.keys():
            t.write("%s\n" %taxa)


def prep_kmfiles(sqzm,taxlist, binlist):
    print("preparing files for keggmapper web tool")
    tax_dict = {}
    bin_dict = {}
    if taxlist:
        tax_list = open(taxlist, "r")
        orftbl = orftable(sqzm)
        for line in tax_list:
            orf = open(orftbl, "r")
            line = line.strip()
            tax_dict[line] = tax(line)
            tax_dict[line].get_ko_list(orf)
            orf.close()
        tax_list.close()
        for taxa in tax_dict.keys():
            tax_dict[taxa].print_ko()

    if binlist:
        bin_list = open(binlist, "r")
        orftbl = orftable(sqzm)
        bintbl = bincontig(sqzm)
        for line in bin_list:
            orf = open(orftbl, "r")
            bin = open(bintbl, "r")
            line = line.strip()
            bin_dict[line] = mag(line)
            bin_dict[line].get_ko_list(bin, orf)
            orf.close()
            bin.close()
        bin_list.close()
        for bin in bin_dict.keys():
            bin_dict[bin].print_ko()
