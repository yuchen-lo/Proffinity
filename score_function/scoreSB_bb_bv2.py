import argparse
import sys
from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt, log, exp
import numpy as np
import pandas as pd
import math

#################
e={}
e[("a","Backbone")]=-0.52
e[("a","ALA")]=-0.19
e[("a","ARG")]=-0.11
e[("a","ASN")]=0.06
e[("a","ASP")]=-0.08
e[("a","CYS")]=0.27
e[("a","GLN")]=-0.12
e[("a","GLU")]=-0.22
e[("a","GLY")]=1.49
e[("a","HIS")]=0.01
e[("a","ILE")]=0.18
e[("a","LEU")]=-0.08
e[("a","LYS")]=-0.15
e[("a","MET")]=-0.09
e[("a","PHE")]=0.09
e[("a","PRO")]=0.75
e[("a","SER")]=0.03
e[("a","THR")]=0.15
e[("a","TRP")]=0.02
e[("a","TYR")]=0.21
e[("a","VAL")]=0.33

e[("b","Backbone")]=0.01
e[("b","ALA")]=0.4
e[("b","ARG")]=0.03
e[("b","ASN")]=0.12
e[("b","ASP")]=0.27
e[("b","CYS")]=-0.34
e[("b","GLN")]=0.07
e[("b","GLU")]=0.41
e[("b","GLY")]=1.45
e[("b","HIS")]=-0.09
e[("b","ILE")]=-0.33
e[("b","LEU")]=0.03
e[("b","LYS")]=0.17
e[("b","MET")]=-0.03
e[("b","PHE")]=-0.26
e[("b","PRO")]=3.41
e[("b","SER")]=-0.09
e[("b","THR")]=-0.27
e[("b","TRP")]=-0.13
e[("b","TYR")]=-0.33
e[("b","VAL")]=-0.4

e[("bp","Backbone")]=0.69
e[("bp","ALA")]=0.06
e[("bp","ARG")]=0.21
e[("bp","ASN")]=0.19
e[("bp","ASP")]=-0.14
e[("bp","CYS")]=0.1
e[("bp","GLN")]=0.21
e[("bp","GLU")]=0.22
e[("bp","GLY")]=1.36
e[("bp","HIS")]=0.18
e[("bp","ILE")]=0.5
e[("bp","LEU")]=0.04
e[("bp","LYS")]=0.16
e[("bp","MET")]=0.35
e[("bp","PHE")]=0.28
e[("bp","PRO")]=-0.71
e[("bp","SER")]=-0.07
e[("bp","THR")]=0.15
e[("bp","TRP")]=0.09
e[("bp","TYR")]=0.19
e[("bp","VAL")]=0.45

e[("al","Backbone")]=1.72
e[("al","ALA")]=1.19
e[("al","ARG")]=0.75
e[("al","ASN")]=-0.78
e[("al","ASP")]=-0.01
e[("al","CYS")]=0.7
e[("al","GLN")]=0.48
e[("al","GLU")]=0.93
e[("al","GLY")]=-1.1
e[("al","HIS")]=-0.07
e[("al","ILE")]=2.58
e[("al","LEU")]=1.39
e[("al","LYS")]=0.31
e[("al","MET")]=1.15
e[("al","PHE")]=1.29
e[("al","PRO")]=4.58
e[("al","SER")]=0.69
e[("al","THR")]=2.07
e[("al","TRP")]=1.28
e[("al","TYR")]=0.85
e[("al","VAL")]=2.76

e[("bl","Backbone")]=2.55
e[("bl","ALA")]=1.21
e[("bl","ARG")]=1.47
e[("bl","ASN")]=0.94
e[("bl","ASP")]=0.97
e[("bl","CYS")]=1.14
e[("bl","GLN")]=1.52
e[("bl","GLU")]=1.44
e[("bl","GLY")]=-1.45
e[("bl","HIS")]=1.49
e[("bl","ILE")]=3.11
e[("bl","LEU")]=2.1
e[("bl","LYS")]=1.34
e[("bl","MET")]=3.67
e[("bl","PHE")]=2.08
e[("bl","PRO")]=4.22
e[("bl","SER")]=0.61
e[("bl","THR")]=1.27
e[("bl","TRP")]=2.36
e[("bl","TYR")]=1.83
e[("bl","VAL")]=2.47

e[("aaa","Backbone")]=-0.82
e[("aaa","ALA")]=-0.26
e[("aaa","ARG")]=-0.12
e[("aaa","ASN")]=0.33
e[("aaa","ASP")]=0.25
e[("aaa","CYS")]=0.33
e[("aaa","GLN")]=-0.13
e[("aaa","GLU")]=-0.24
e[("aaa","GLY")]=1.51
e[("aaa","HIS")]=0.04
e[("aaa","ILE")]=0.1
e[("aaa","LEU")]=-0.15
e[("aaa","LYS")]=-0.09
e[("aaa","MET")]=-0.15
e[("aaa","PHE")]=0.06
e[("aaa","PRO")]=2.32
e[("aaa","SER")]=0.27
e[("aaa","THR")]=0.35
e[("aaa","TRP")]=-0.01
e[("aaa","TYR")]=0.21
e[("aaa","VAL")]=0.28

e[("bbb","Backbone")]=-0.57
e[("bbb","ALA")]=0.31
e[("bbb","ARG")]=0.39
e[("bbb","ASN")]=0.81
e[("bbb","ASP")]=1.19
e[("bbb","CYS")]=-0.31
e[("bbb","GLN")]=0.43
e[("bbb","GLU")]=0.57
e[("bbb","GLY")]=1.44
e[("bbb","HIS")]=0.28
e[("bbb","ILE")]=-0.49
e[("bbb","LEU")]=-0.16
e[("bbb","LYS")]=0.58
e[("bbb","MET")]=-0.13
e[("bbb","PHE")]=-0.37
e[("bbb","PRO")]=4.44
e[("bbb","SER")]=0.12
e[("bbb","THR")]=-0.16
e[("bbb","TRP")]=-0.28
e[("bbb","TYR")]=-0.45
e[("bbb","VAL")]=-0.57

e[("bpbbp","Backbone")]=-0.18
e[("bpbbp","ALA")]=0.5
e[("bpbbp","ARG")]=-0.59
e[("bpbbp","ASN")]=0.8
e[("bpbbp","ASP")]=1.44
e[("bpbbp","CYS")]=-0.42
e[("bpbbp","GLN")]=-0.3
e[("bpbbp","GLU")]=0.13
e[("bpbbp","GLY")]=2.07
e[("bpbbp","HIS")]=0.22
e[("bpbbp","ILE")]=-0.31
e[("bpbbp","LEU")]=0.01
e[("bpbbp","LYS")]=-0.39
e[("bpbbp","MET")]=-0.42
e[("bpbbp","PHE")]=0.25
e[("bpbbp","PRO")]=3.4
e[("bpbbp","SER")]=0.33
e[("bpbbp","THR")]=-0.01
e[("bpbbp","TRP")]=0.36
e[("bpbbp","TYR")]=0.52
e[("bpbbp","VAL")]=-0.56

e[("bbbp","Backbone")]=-0.11
e[("bbbp","ALA")]=0.46
e[("bbbp","ARG")]=-0.01
e[("bbbp","ASN")]=0.71
e[("bbbp","ASP")]=1.35
e[("bbbp","CYS")]=-0.11
e[("bbbp","GLN")]=-0.05
e[("bbbp","GLU")]=0.28
e[("bbbp","GLY")]=1.75
e[("bbbp","HIS")]=0.13
e[("bbbp","ILE")]=-0.48
e[("bbbp","LEU")]=-0.08
e[("bbbp","LYS")]=0.04
e[("bbbp","MET")]=-0.16
e[("bbbp","PHE")]=-0.4
e[("bbbp","PRO")]=3.32
e[("bbbp","SER")]=0.4
e[("bbbp","THR")]=-0.1
e[("bbbp","TRP")]=-0.08
e[("bbbp","TYR")]=-0.38
e[("bbbp","VAL")]=-0.44

e[("bpbb","Backbone")]=-0.04
e[("bpbb","ALA")]=0.52
e[("bpbb","ARG")]=-0.13
e[("bpbb","ASN")]=0.75
e[("bpbb","ASP")]=1.47
e[("bpbb","CYS")]=-0.42
e[("bpbb","GLN")]=-0.21
e[("bpbb","GLU")]=0.91
e[("bpbb","GLY")]=1.81
e[("bpbb","HIS")]=0.06
e[("bpbb","ILE")]=-0.53
e[("bpbb","LEU")]=-0.02
e[("bpbb","LYS")]=0.08
e[("bpbb","MET")]=-0.15
e[("bpbb","PHE")]=-0.1
e[("bpbb","PRO")]=3.16
e[("bpbb","SER")]=0.34
e[("bpbb","THR")]=-0.29
e[("bpbb","TRP")]=-0.15
e[("bpbb","TYR")]=-0.17
e[("bpbb","VAL")]=-0.63

e[("bpbpbp","Backbone")]=0.44
e[("bpbpbp","ALA")]=-0.04
e[("bpbpbp","ARG")]=0.28
e[("bpbpbp","ASN")]=0.86
e[("bpbpbp","ASP")]=0.6
e[("bpbpbp","CYS")]=0.05
e[("bpbpbp","GLN")]=0.01
e[("bpbpbp","GLU")]=0.25
e[("bpbpbp","GLY")]=1.92
e[("bpbpbp","HIS")]=0.32
e[("bpbpbp","ILE")]=0.16
e[("bpbpbp","LEU")]=0.05
e[("bpbpbp","LYS")]=0.06
e[("bpbpbp","MET")]=0.25
e[("bpbpbp","PHE")]=0.93
e[("bpbpbp","PRO")]=-0.82
e[("bpbpbp","SER")]=0.77
e[("bpbpbp","THR")]=0.77
e[("bpbpbp","TRP")]=2.39
e[("bpbpbp","TYR")]=0.78
e[("bpbpbp","VAL")]=0.55


def yieldPhiPsi(model):
    for chain in model:
        polypeptides = PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly):
                phi, psi = phi_psi[res_index]
                yield dict(chain_id=chain.id, 
                           poly_index=poly_index, 
                           residue_name=residue.resname, 
                           author_residue_number=residue.id[1], 
                           phi_radian=phi, psi_radian=psi)


def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, then it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle


def conformational_state (phi, psi) :
    state="NaN"
    if phi <= 0 and psi >= -120 and psi <= 60:
        state="a"
        return state
    elif phi < -90 and state != "a":
        state="b"
        return state
    elif phi >= -90 and phi <= 0 and state != "a":
        state="bp"
        return state
    elif phi > 0 and psi >= -60 and psi < 120:
        state="al"
        return state
    elif phi > 0 and state != "al":
        state="bl"
        return state
    else:
        return state
    
def triple_conformational_state (l,m,r) :
    t_state="NaN"
    if l == "b" and m == "b" and r == "b":
        t_state=l+m+r
        return t_state
    elif l == "bp" and m == "b" and r == "bp":
        t_state=l+m+r
        return t_state
    elif l == "b" and m == "b" and r == "bp":
        t_state=l+m+r
        return t_state
    elif l == "bp" and m == "b" and r == "b":
        t_state=l+m+r
        return t_state
    elif l == "bp" and m == "bp" and r == "bp":
        t_state=l+m+r
        return t_state
    elif l == "a" and m == "a" and r == "a":
        t_state=l+m+r
        return t_state
    else:
        return t_state
    
def scoreSB_BB_B(pdbid, ctype):

    p = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    s = p.get_structure("X", "./model/"+pdbid+".pdb")
    #s = p.get_structure("X", "./modeller_output_ab/best_model_"+args.pdbid+".pdb")

    #generate phi psi angles
    df = pd.DataFrame(yieldPhiPsi(s[0]))

    df['phi_degree'] = df.phi_radian.apply(degrees)
    df['psi_degree'] = df.psi_radian.apply(degrees)

    #assign secondary conformation
    df['state']="NaN"
    df['t_state']="NaN"

    pd.options.mode.chained_assignment = None

    for index, row in df.iterrows():
        #print(row['phi_degree'], row['psi_degree'])
        if np.isnan(row['phi_degree']) | np.isnan(row['psi_degree']):
            continue
        else:
            #print("conformational state:"+str(conformational_state(row['phi_degree'],row['psi_degree'])))
            df['state'][index]=conformational_state(row['phi_degree'],row['psi_degree'])

    for index, row in df.iterrows():
        l=index
        m=index+1
        r=index+2
        
        if r <= df.shape[0] and df['state'][m] != "NaN":
            #print("r="+str(r))
            df['t_state'][m]=triple_conformational_state(df['state'][l],df['state'][m],df['state'][r])

    #compute total secondary energy
    df['t_energy']=0
    df['energy']=0

    for index, row in df.iterrows():
        l=index
        m=index+1
        r=index+2

        if r <= df.shape[0] and df['state'][m] != "NaN":
            #print("state="+df['state'][m])
            df['energy'][m]= e[(df['state'][m],df['residue_name'][m])]
            
    
        
        if r <= df.shape[0] and df['t_state'][m] != "NaN":
            df['t_energy'][m]= e[(df['t_state'][m],"Backbone")] + e[(df['t_state'][m],df['residue_name'][m])]+ e[(df['t_state'][m],df['residue_name'][l])]+ e[(df['t_state'][m],df['residue_name'][r])]
        

    #extract energy from the interface
    file1 = open('./raw_graphv2/'+pdbid+'_CA_output_bb_b.txt', 'r')

    Lines = file1.readlines()
    df2=pd.DataFrame(Lines)
    df2.columns=['col']
    idx=list(df2['col'].str.split(',').str[2].unique().astype(int))

    with open('./raw_graphv2/' + pdbid + '_CA_output_bb_b.txt', 'w') as f:
        print("")

    for line in Lines:
        chain_A=line.split(',')[0]
        resA=line.split(',')[1]
        resAid=line.split(',')[2]
        chain_B=line.split(',')[3]
        resB=line.split(',')[4]
        dist=line.split(',')[6]

        #print("chain A:"+chain_A+" chain_B:"+chain_B)
        if (chain_A==chain_B):
            t_eng_i=df['t_energy'][(df['author_residue_number']==int(resAid))].values[0]
            eng_i=df['energy'][(df['author_residue_number']==int(resAid))].values[0]
            line=line.strip("\n")+","+str(round(eng_i,2))+","+str(round(t_eng_i,2))
            with open('./raw_graphv2/' + pdbid + '_CA_output_bb_b.txt', 'a') as f:
                print(line, file=f)


    t_eng=df[df['author_residue_number'].isin(idx)]['t_energy'].sum()
    eng=df[df['author_residue_number'].isin(idx)]['energy'].sum()


    print("scoreSB_BB_tB:",str(round(t_eng,2))," ctype:",ctype)
    print("scoreSB_BB_B:",str(round(eng,2))," ctype:",ctype)
    return round(t_eng,2),round(eng,2)




