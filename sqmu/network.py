from collections import defaultdict
from collections import Counter
import math
import numpy as np
from sqmu.sqmfiles import orftable, projname, sqmsamples


class node_ko(object):
    def __init__(self, ko_id, samples):
        self.ko_id = ko_id
        self.kofun = "n.d."
        self.koclass = "n.d."
        self.cogfun = "n.d."
        self.cogclass = "n.d."
        self.tpm = np.zeros(len(samples))

    def add_info(self, annotation):
        tmparr = []
        for i in range(0, len(self.tpm)):
            a = annotation[16+i]
            if len(a) > 0:
                tmparr.append(a)
            else:
                tmparr.append("0")
        self.tpm = self.tpm + [float(a) for a in tmparr]

        if "n.d." in self.cogfun and len(annotation[13]) > 0:
            self.cogfun = annotation[13]
        if "n.d." in self.kofun and len(annotation[10]) > 0:
            self.kofun = annotation[10]
        if "n.d." in self.koclass and len(annotation[11]) > 0:
            self.koclass = annotation[11]
        if "n.d." in self.cogclass and len(annotation[14]) > 0:
            self.cogclass = annotation[14]


class node_org(object):
    def __init__(self, org_id, samples):
        self.org_id = org_id
        self.orgset = set()
        self.tpm = np.zeros(len(samples))

    def add_info(self, annotation):
        tmparr = []
        for i in range(0, len(self.tpm)):
            a = annotation[16+i]
            if len(a) > 0:
                tmparr.append(a)
            else:
                tmparr.append("0")
        self.tpm = self.tpm+[float(a) for a in tmparr]
        self.orgset.add(annotation[8])


class edge_interaction(object):
    def __init__(self, edge_id, samples):
        self.edge_id = edge_id
        self.tpm = np.zeros(len(samples))

    def add_info(self, annotation):
        tmparr = []
        for i in range(0, len(self.tpm)):
            a = annotation[16+i]
            if len(a) > 0:
                tmparr.append(a)
            else:
                tmparr.append("0")
        self.tpm = self.tpm + [float(a) for a in tmparr]


def printnetwork(ko_dict, org_dict, edge_dict, samples, out):
    # Create output files node annotation(ko and org)
    file_node_ko = open(out+"_ko.tsv", "w")
    file_node_org = open(out+"_org.tsv", "w")
    file_node_ko.write("node_id\tabundance\tcog_func\tcog_class\tko_func\tko_class\ttype\n")
    file_node_org.write("node_id\tabundance\tall_organisms\ttype\n")
    for ko in ko_dict.keys():
        # print ko annotation table
        file_node_ko.write("%s\t%s\t%s\t%s\t%s\t%s\tko\n"
                           % (ko, ko_dict[ko].tpm[0], ko_dict[ko].cogfun, ko_dict[ko].cogclass,
                              ko_dict[ko].kofun, ko_dict[ko].koclass))
    for org in org_dict.keys():
        # print org annotation table
        file_node_org.write("%s\t%s\t[%s]\torg\n" % (org, org_dict[org].tpm[0],
                                                     ",".join(org_dict[org].orgset)))
    file_node_ko.close()
    file_node_org.close()
    for i,sample in enumerate(samples):
        with open(out+"_"+sample + "_edge.tsv", "w") as file_edge:
            file_edge.write("source\tinteraction\ttarget\tabundance\tlog2_abundance\n")
            for edge in edge_dict.keys():
                ko, org = edge.split("\t")
                # print edge info for each sample
                if edge_dict[edge].tpm[i] > 0:
                    log2 = math.log(edge_dict[edge].tpm[i]+2, 2)
                    file_edge.write("%s\t%s\t%s\t%s\t%s\n" %
                                    (ko, sample, org, edge_dict[edge].tpm[i], log2))


def build_network(sqm, kolist, taxlevel, out):
    print("buildingn network taxonomic level")
    ko_dict = {}
    org_dict = {}
    edge_dict = {}
    sqm_orf = orftable(sqm)
    samples = sqmsamples(sqm)
    ko_set = set()
    abundance = defaultdict(list)
    # Read file 1 with ko list and save into ko_set variable
    with open(kolist, "r") as ko_file:
        for line in ko_file:
            ko = line.rstrip()
            ko_set.add(ko)
    # Read file 13 orftable build network representation
    with open(sqm_orf, "r") as orf_file:
        next(orf_file)
        first_line = next(orf_file)  # header is now in second line v1.0
        for line in orf_file:
            orf_annot = line.split("\t")
            # for v1 colnum 9-tax 10-kegid 11-keggfun 12-kegclass 15-COGclass
            # [17, 18,...] TPM [sample1, sample2,...]
            # index=colnum-1
            if len(orf_annot[9]) > 1 and orf_annot[9].replace("*", "") in ko_set:  # ko in list
                ko = orf_annot[9]
                ko = ko.replace("*", "")  # replace *
                taxonomy = orf_annot[8].split(";")
                if taxlevel not in orf_annot[8]:
                    if len(orf_annot[8]) == 0:
                        org = "Unclassified"
                    else:
                        org = taxonomy[-1]  # select last taxonomy level
                else:
                    for level in taxonomy:  # select specified taxonomy level
                        if taxlevel in level:
                            org = level
                if org not in org_dict.keys():
                    org_dict[org] = node_org(org, samples)
                org_dict[org].add_info(orf_annot)
                if ko not in ko_dict.keys():
                    ko_dict[ko] = node_ko(ko, samples)
                ko_dict[ko].add_info(orf_annot)
                edge_id = ko+"\t"+org
                if edge_id not in edge_dict.keys():
                    edge_dict[edge_id] = edge_interaction(edge_id, samples)
                edge_dict[edge_id].add_info(orf_annot)
    printnetwork(ko_dict, org_dict, edge_dict, samples, out)
