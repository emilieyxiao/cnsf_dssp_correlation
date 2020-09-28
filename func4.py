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
    #Reading File
    source = open("/home/emiliex/projects/ctb-yxia/emiliex/cnsf/"+ consf_file, "r")
    content = source.read()
    source.close()
    #Stripping lines at beginning and end
    list = content.split('\n')
    list = list[15:]
    list = list[:-5]
    #Split makes sublists (each row in list becomes a separate list)
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
    #print('Cnsf length:')
    #print(len(evoln))
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
    #Main loop
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

"""
a,b,l = extractdssp("1nex.dssp","B")
print(b)
print(l)
"""

def corrfunc(cf,df):
    import numpy as np
    #Solving chain's correlation
    evo = []
    acc = []
    for x in cf:
        if x in df:
            acc.append(df[x])
            evo.append(cf[x])
    #print(acc)
    #print(evo)
    #print(len(acc))
    #print(len(evo))
    co = np.corrcoef(evo,acc)
    sc = co[0,1]
    return sc

"""
    #Finding chain break
    count = -1
    chain_break = []
    for line in lines:
        count = count+1            #First line is 0 and so on...
        if '!*' in line:
            chain_break.append(count)
    chain_break.append(len(lines)) #Last Line
    

    #Storing acc scores in separate dictionaries for each chain
    dic = {}
    start = start_line
    for i in range(0,len(chain_break)):
        end = chain_break[i]
        chain = lines[start][11]
        if chain == desired_chain:
            for line in lines[start:end]:              
                r = line[6:10].strip()          # position is 6th-10th character
                if r in cnsf:
                    score = line[35:38].strip()         # score is 6th-10th character
                    aa = norm.get(line[13].strip())     # get norm score from legend based on amino acid
                    #Normalizing
                    if isinstance(n,int):
                        dic[r] = (float(score)/aa)
                    else:                                                                                                              dic[r] = float(score)/6         #For cysteine bridge
        #Iterating for a next chain
        start = end
        i = i+1
    print('Solvent accessibility:')
    print(dic)
    return dic
"""

