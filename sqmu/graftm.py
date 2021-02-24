'''Read graftm results
find orf names in results, count tpm for each detected
orf with homology, add to abundance and specified taxonomy
'''
from collections import defaultdict
from sqmu.sqmfiles import projname, orftable, sqmsamples
import numpy as np
import os
import re

def graftm_samplename(graft):
    graftpath=os.path.realpath(graft)
    combinedct=os.path.join(graftpath,"combined_count_table.txt")
    with open(combinedct,"r") as f:
        header=next(f)
    samplenames=header.split("\t")
    print(samplenames[1:-1])
    return(samplenames[1:-1])

def graftm_readtax(graft, sample):
    readtax=defaultdict(str)
    filename=sample+"_read_tax.tsv"
    graftpath=os.path.realpath(graft)
    readtaxpath=os.path.join(graftpath,sample)
    filename=os.path.join(readtaxpath,filename)
    with open(filename, "r") as f:
        for line in f:
            line=line.rstrip()
            read,tax=line.split("\t")[0:2]
            read="_".join(read.split("_")[0:-3])
            readtax[read]=tax
    return(readtax)

def graftm_abundance(sqm,graft,ab,out):
    ab=16
    print("Extracting abundance values for graftM results")
    ab_counttable={}
    sqm_orf=orftable(sqm)
    sqm_samples=sqmsamples(sqm)
    graftsamples=graftm_samplename(graft)
    for gsample in graftsamples:
        read_tax_dict=graftm_readtax(graft,gsample)
        with open(sqm_orf, "r") as f:
            next(f)
            header=next(f)
            header=header.split("\t")
            for line in f:
                line=line.split("\t")
                orfid=line[0]
                line_ab=[]
                if orfid in read_tax_dict.keys():
                    for i,sqm_sam in enumerate(sqm_samples):
                        value=float(line[ab+i])
                        line_ab.append(value)
                    line_ab=np.array(line_ab)
                    tax=read_tax_dict[orfid]
                    if tax not in ab_counttable.keys():
                        ab_counttable[tax]=np.zeros(len(sqm_samples))
                        ab_counttable[tax]=ab_counttable[tax]+line_ab
                    else:
                        ab_counttable[tax]=ab_counttable[tax]+line_ab
        #printcounttable(ab_counttable, sqm_samples, out)
        outfile=out+"_"+gsample+"_combined_count_table.txt"
        with open(outfile, "w") as out:
            out.write("#ID\t")
            for sqmsQam in sqm_samples:
                out.write("%s\t" %sqmsam)
            out.write("ConsensusLineage\n")
            for i,tax in enumerate(sorted(ab_counttable.keys())):
                out.write("%s\t" %i)
                for j in ab_counttable[tax]:
                    out.write("%s\t" %j)
                out.write("%s\n" %tax)
