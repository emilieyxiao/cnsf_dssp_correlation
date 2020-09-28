# functions to extract and normalize dssp and consurf scores then calculate correlation
# use with main.py

def normalize(score,aminoa):
    #Normalization array for residue size
    norm = {'A':115,'R':11,'N':8,'D':8,'C':6,'Q':9,'E':9,'G':91,'H':10,'I':8,'L':159,'K':206,'M':8,'F':11,'P':145,'S':128,'T':7,'W':14,'Y':12,'V':135}
    aa = norm.get(aminoa)
    if isinstance(aa,int):
        nscore = (float(score)/aa)
    else:
        nscore = float(score)/6         #For cysteine bridge
    #print(type(nscore))
    nscore = round(nscore,3)
    #print(type(nscore))
    #nscore = format(nscore,'.3f')
    return nscore

def extractcnsf(consf_file):
    source = open("/home/emiliex/projects/ctb-yxia/emiliex/cnsf/"+ consf_file, "r")
    content = source.read()
    source.close()
    #Stripping lines at beginning and end
    list = content.split('\n')
    list = list[15:]
    list = list[:-5]
    split = [i.split() for i in list]
    evoln = {}
    i = 0
    for i in range(len(split)):
        if ':' in split[i][2]:
            res = split[i][2]
            res = res[3:]
            res = res[:-2]
            score = float(split[i][3])
            evoln[res] = score
            i = i+1
    return evoln

def extractdssp(dssp_file,desired_chain):
    #Reading file
    dssp = open("/home/emiliex/projects/ctb-yxia/emiliex/dssp/"+ dssp_file,'r')
    lines = dssp.readlines()
    dssp.close
    #Finding start line
    count = -1
    for line in lines:
        count = count+1            #First line is 0 and so on...
        if '#  RESIDUE AA STRUCTURE' in line:
            start_line = count + 1
    asa ={}
    z = 0 #counting residues with score = zero
    for line in lines[start_line:]:
        chain = line[11]
        if chain == desired_chain:   #match - right chain
            res = line[6:10].strip()
            sc = line[35:38].strip() #score
            aa = line[13].strip()    #amino acid
            asa[res] = normalize(sc,aa) #store normalized score 
            if float(sc) == 0:
                z = z+1
    #print('DSSP length:')
    #print(len(asa))
    pz = 100*z/len(asa) #percent of zero
    pz = round(pz,3)
    return asa,pz

def corrfunc(cf,df):
    import numpy as np
    #Solving chain's correlation
    evo = []
    acc = []
    for x in cf:
        if x in df:
            acc.append(df[x])
            evo.append(cf[x])
    co = np.corrcoef(evo,acc)
    sc = co[0,1]
    return sc

