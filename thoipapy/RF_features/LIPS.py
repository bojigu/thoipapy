# here is the python script to calculate LIPS score for TMPs, just a convertion from lips.pl

import math

sequence = ""
Lips_input_file = r"/scratch2/zeng/Exp_Pred_161020/a3m/Q9Y286.MEM.lips.input"  ##give the imput alignment without gaps
file = open(Lips_input_file, "r")
sequence = ' '.join(line.strip() for line in file)

n = 0
sump = 0
sumlip = 0
sume = {}  # the sum of entropy for each one of the seven surfaces
sumf = 0
sumim = {}  # the sum of lipophilicity for each surfaces
aanum = {}  # the number of residues for each one of seven surfaces
resnum = 1
amino = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

###propi and propm are both from TMLIP scale, propi means the residue lipophilicity propensity in membrane headgroup region
##and propm is the residue lipophilicity propensity in hydrophobic core region
##in TMLIP scale paper, the membrane headgroup regions is defined as the first (tmlen/5) and (tmlen-tmlen/5) residues which
##more likely the membrane bilayer
##while the other residues are defined as hydrophobic core region


propi = {
    'A': 0.71,
    'R': 1.47,
    'N': 0.96,
    'D': 1.20,
    'C': 1.16,
    'Q': 0.61,
    'E': 0.90,
    'G': 0.48,
    'H': 0.82,
    'I': 1.11,
    'L': 1.18,
    'K': 2.38,
    'M': 1.38,
    'F': 1.57,
    'P': 0.99,
    'S': 0.69,
    'T': 0.72,
    'W': 2.45,
    'Y': 1.23,
    'V': 0.98
}

propm = dict(
    A=0.82,
    R=0.18,
    N=0.19,
    D=0.29,
    C=1.01,
    Q=0.26,
    E=0.19,
    G=0.35,
    H=0.12,
    I=1.88,
    L=1.71,
    K=0.42,
    M=1.02,
    F=1.97,
    P=0.65,
    S=0.55,
    T=0.66,
    W=1.65,
    Y=0.94,
    V=1.77
)

tmp = sequence.split(' ')
nrow = len(tmp)
ncol = len(tmp[0])
bnum = ncol / 5
oc = {}
prob = {}
entropy = {}  ##residue entropy
exp_entropy = {}  # the final exponential entropy
lips = {}  ##the lipophilicity score

for i in range(nrow):
    for j in range(ncol):
        residue = tmp[i][j]
        res_j = ' '.join((residue, str(j)))
        if (res_j in oc.keys()):
            oc[res_j] = oc[res_j] + 1
        else:
            oc[res_j] = 1

for j in range(ncol):
    for res in amino:
        if (' '.join((res, str(j))) in oc):
            prob[res] = oc[' '.join((res, str(j)))] / nrow
            if (j in entropy.keys()):
                entropy[j] = entropy[j] + prob[res] * math.log(prob[res])  # the entropy calculation
            else:
                entropy[j] = prob[res] * math.log(prob[res])
            if ((j <= bnum) or (j > ncol - bnum)):  ###here is the membrane headgroup residues
                if (j in lips.keys()):
                    lips[j] = lips[j] + prob[res] * propi[res]
                else:
                    lips[j] = prob[res] * propi[res]
            else:  ###here is the hydrophobic region residues
                if (j in lips.keys()):
                    lips[j] = lips[j] + prob[res] * propm[res]
                else:
                    lips[j] = prob[res] * propm[res]
    exp_entropy[j] = 2.718 ** ((-1) * entropy[j])  # expontional entropy

for j in sorted(exp_entropy):
    res = tmp[0][j]
    m = resnum + j
    sump = sump + exp_entropy[j]
    sumlip = sumlip + lips[j]

for i in range(4):  # for the first 4 surfaces
    print("SURFACE", "%s" % i, ":")
    j = i
    while j < ncol:
        res = tmp[0][j]
        if (i in sumim.keys()):
            sumim[i] = sumim[i] + lips[j]  # the sum of lipophilicity for surface i
        else:
            sumim[i] = lips[j]
        prop = lips[j]
        if (i in sume.keys()):
            sume[i] = sume[i] + exp_entropy[j]  # the sum of entropy for surface i
        else:
            sume[i] = exp_entropy[j]
        if (i in aanum.keys()):
            aanum[i] = aanum[i] + 1  # the sum of the residue numbers for surface i
        else:
            aanum[i] = 1
        rn = j + resnum
        # r3=residuename123(res)
        print("%3s" % rn, res, "%6.3f" % prop,
              "%6.3f" % exp_entropy[j])  # print residue information which is in surface i
        k = j + 3
        while (k <= j + 4):  # here add the the residues of i+3 and i+4 into surface i to form heptad repeat
            if (k < ncol):
                res = tmp[0][k]
                # r3=residuename123(res)
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[k]
                else:
                    sumim[i] = lips[k]
                prob = lips[k]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[k]
                else:
                    sume[i] = exp_entropy[k]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1
                else:
                    aanum[i] = 1
                rn = k + resnum
                print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
            k = k + 1
        j = j + 7
for i in range(4, 7):  # for surfaces from 4 to 6
    print("SURFACE", "%s" % i, ":")
    j = i
    while j < ncol:
        res = tmp[0][j]
        if (i in sumim.keys()):
            sumim[i] = sumim[i] + lips[j]
        else:
            sumim[i] = lips[j]
        prob = lips[j]
        if (i in sume.keys()):
            sume[i] = sume[i] + exp_entropy[j]
        else:
            sume[i] = exp_entropy[j]
        if (i in aanum.keys()):
            aanum[i] = aanum[i] + 1
        else:
            aanum[i] = 1
        rn = j + resnum
        # r3=residuename123(res)
        print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[j])
        k = j + 3
        while (k <= j + 4):
            if (k < ncol):
                res = tmp[0][k]
                # r3=residuename123(res)
                if (i in sumim.keys()):
                    sumim[i] = sumim[i] + lips[k]
                else:
                    sumim[i] = lips[k]
                prob = lips[k]
                if (i in sume.keys()):
                    sume[i] = sume[i] + exp_entropy[k]
                else:
                    sume[i] = exp_entropy[k]
                if (i in aanum.keys()):
                    aanum[i] = aanum[i] + 1
                else:
                    aanum[i] = 1
                rn = k + resnum
                print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
            k = k + 1
        j = j + 7
    k = i - 4
    while (k <= i - 3):  # here adding residues at the first 7 positions
        if (k < ncol):
            res = tmp[0][k]
            # r3=residuename123(res)
            if (i in sumim.keys()):
                sumim[i] = sumim[i] + lips[k]
            else:
                sumim[i] = lips[k]
            prob = lips[k]
            if (i in sume.keys()):
                sume[i] = sume[i] + exp_entropy[k]
            else:
                sume[i] = exp_entropy[k]
            if (i in aanum.keys()):
                aanum[i] = aanum[i] + 1
            else:
                aanum[i] = 1
            rn = k + resnum
            print("%3s" % rn, res, "%6.3f" % prob, "%6.3f" % exp_entropy[k])
        k = k + 1
print("SURFACE LIPOPHILICITY ENTROPY   LIPS")
for i in sumim.keys():
    avpim = sumim[i] / aanum[i]  # average lipophilicity for surface i
    avpim = avpim * 2
    ave = sume[i] / aanum[i]  # average entropy for surface i
    peim = avpim * ave  # average entropy*lipophilicity for surface i which is LIPS score
    print("%s" % i, "%10.3f" % avpim, "%8.3f" % ave,
          "%8.3f" % peim)  # print seven surfaces and see which surface with lowewst LIPS score
