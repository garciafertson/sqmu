from collections import defaultdict

def run_kmod_kmap(km_files,modlist,outname):
    print("buil matrix from keggmaper result files")
    module_dict={}
    module_ls=[]
    outname=outname+"_keggmap.tsv"
    table_dict=defaultdict(list)
    tax_ls=[]
    categories_set={"complete)":3, "2 blocks missing)":1, "incomplete)":0, "1 block missing)":2}

    with open(modlist) as f:
        for line in f:
            line=line.strip()
            mod=line.split("\t")[0]
            module_dict[mod]=line.split("\t")[1]
            module_ls.append(mod)
    #print(module_ls)
    for kmresult in km_files:
        if len(kmresult.split(".")) > 0: tax=".".join(kmresult.split(".")[0:-1])
        else: tax=kmresult
        tax_ls.append(tax)
        with open(kmresult, "r") as kmap:
            for line in kmap:
                line=line.rstrip()
                for module in module_ls:
                    if module in line:
                        km=line.split("(")
                        table_dict[module,tax].append(categories_set[km[-1]])

    with open(outname,"w") as f:
        f.write("taxid")
        for module in module_ls:
            print(module_dict[module])
            f.write("\t%s" %module_dict[module])
        f.write("\n")
        for tax in tax_ls:
            f.write("%s" %tax)
            for module in module_ls:
                f.write("\t%s" %table_dict[module,tax][0])
            f.write("\n")
