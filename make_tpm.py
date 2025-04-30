import os
import sys

#run with python make_rna_ss.py path/to/directory/with/kallisto_outputs

mydir = sys.argv[1]
if not os.path.exists(mydir):
    sys.exit("{0} is not a valid directory".format(mydir))

gene_exp_dict = {}
sample_list = []
mysamples = os.listdir(mydir)
for asample in mysamples:
    if not os.path.isdir(mydir + "/" + asample): continue
    if not os.path.exists(mydir + "/" + asample + "/abundance.tsv"): continue
    fh = open(mydir + "/" + asample + "/abundance.tsv")
    sample_list.append(asample)
    fh.readline()
    for x in fh:
        y = x.strip().split('\t')
        mygene = y[0]
        mytpm = float(y[-1])
        if not mygene in gene_exp_dict: gene_exp_dict[mygene] = {}
        gene_exp_dict[mygene][asample] = mytpm
    fh.close()
fh = open("merged_gene_tpms.csv2",'w')
myheader = ["GeneID"] + sorted(sample_list)
#print(len(gene_exp_dict))
#print(list(gene_exp_dict)[:10])
fh.write(",".join(myheader)+"\n")
for agene in sorted(list(gene_exp_dict)):
    plist = [agene]
    for asample in sorted(sample_list):
        plist.append(gene_exp_dict[agene][asample])
    fh.write(",".join(map(str,plist))+"\n")
    
