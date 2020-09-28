import os
import numpy as np
import csv
from func4 import extractdssp,extractcnsf,corrfunc
 
#Opening files
dssp = [i for i in os.listdir("/home/emiliex/projects/ctb-yxia/emiliex/dssp") if i.endswith("dssp")]
cnsf = [i for i in os.listdir("/home/emiliex/projects/ctb-yxia/emiliex/cnsf") if i.endswith("consurf_summary.txt")]

c = {}
for file in cnsf:
    c[file] = file[:5],file[:4].lower(),file[4]
    # fullname,residue, chain
d = {}
for file in dssp:
    d[file] = file[:4].lower()
    # residue

residue = {}
score = {}
percent = {}
chai = {}

for i in c:
    name = c[i][0]
    res = c[i][1]
    for j in d:
        if d[j] == res:  #if residues match:
            chain = c[i][2]
            cnsf = extractcnsf(i)
            dssp,pz = extractdssp(j,chain)
            residue[name] = res
            score[name] = corrfunc(cnsf,dssp)
            percent[name] = pz
            chai[name] = chain
            break
print(len(c))
#Output to csv
dicts = residue,score,percent,chai
with open('fullresults.csv','wb') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['Name','Residue','Correlation','Percent 0 Exposure','Chain'])#headers
    for key in residue.iterkeys():
        writer.writerow([key]+[d[key] for d in dicts])


