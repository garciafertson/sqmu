#read sqm files
import re
import os


def projname(sqzm):
    sqmpath=os.path.realpath(sqzm)
    configfile=os.path.join(sqmpath,"SqueezeMeta_conf.pl")
    with open(configfile,"r") as f:
        for line in f:
            if re.search(r'^\$projectname\s=\s\"(?P<name>\w+)\"',line):
                projectname=re.search(r'^\$projectname\s=\s\"(?P<name>\w+)\"',line).group('name')
    return(projectname)

def orftable(sqzm):
    projectname=projname(sqzm)
    sqmpath=os.path.realpath(sqzm)
    orftable="13."+projectname+".orftable"
    orftable=os.path.join("results",orftable)
    orftable=os.path.join(sqmpath,orftable)
    return(orftable)

def bincontig(sqzm):
    projectname=projname(sqzm)
    sqmpath=os.path.realpath(sqzm)
    contigtable="20."+projectname+".contigtable"
    contigtable=os.path.join("results",contigtable)
    contigtable=os.path.join(sqmpath,contigtable)
    return(contigtable)

def sqmsamples(sqzm):
    projectname=projname(sqzm)
    sqmpath=os.path.realpath(sqzm)
    mapstat="10."+projectname+".mappingstat"
    mapstat=os.path.join("results", mapstat)
    mapstat=os.path.join(sqmpath,mapstat)
    samples=set()
    with open(mapstat) as f:
        for line in f:
            if not re.match(r'^#', line):
                samples.add(line.split("\t")[0])
    return(samples)
