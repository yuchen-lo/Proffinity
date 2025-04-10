import argparse
import sys
from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt, log

#################
#Dosztanyi et al. (2005)
#Eij=eij+err-eir-ejr
eij={}
eij[("ALA","ALA")]=-1.65
eij[("ALA","CYS")]=eij[("CYS","ALA")]=-2.83
eij[("ALA","ASP")]=eij[("ASP","ALA")]=1.16
eij[("ALA","GLU")]=eij[("GLU","ALA")]=1.8
eij[("ALA","PHE")]=eij[("PHE","ALA")]=-3.73
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-0.41
eij[("ALA","HIS")]=eij[("HIS","ALA")]=1.9
eij[("ALA","ILE")]=eij[("ILE","ALA")]=-3.69
eij[("ALA","LYS")]=eij[("LYS","ALA")]=0.49
eij[("ALA","LEU")]=eij[("LEU","ALA")]=-3.01
eij[("ALA","MET")]=eij[("MET","ALA")]=-2.08
eij[("ALA","ASN")]=eij[("ASN","ALA")]=0.66
eij[("ALA","PRO")]=eij[("PRO","ALA")]=1.54
eij[("ALA","GLN")]=eij[("GLN","ALA")]=1.20
eij[("ALA","ARG")]=eij[("ARG","ALA")]=0.98
eij[("ALA","SER")]=eij[("SER","ALA")]=-0.08
eij[("ALA","THR")]=eij[("THR","ALA")]=0.46
eij[("ALA","VAL")]=eij[("VAL","ALA")]=-2.31
eij[("ALA","TRP")]=eij[("TRP","ALA")]=0.32
eij[("ALA","TYR")]=eij[("TYR","ALA")]=-4.62
eij[("CYS","CYS")]=-39.58
eij[("CYS","ASP")]=eij[("ASP","CYS")]=-0.82
eij[("CYS","GLU")]=eij[("GLU","CYS")]=-0.53
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-3.07
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-2.96
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-4.98
eij[("CYS","ILE")]=eij[("ILE","CYS")]=0.34
eij[("CYS","LYS")]=eij[("LYS","CYS")]=-1.38
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-2.15
eij[("CYS","MET")]=eij[("MET","CYS")]=1.43
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-4.18
eij[("CYS","PRO")]=eij[("PRO","CYS")]=-2.13
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-2.91
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-0.41
eij[("CYS","SER")]=eij[("SER","CYS")]=-2.33
eij[("CYS","THR")]=eij[("THR","CYS")]=-1.84
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-0.16
eij[("CYS","TRP")]=eij[("TRP","CYS")]=4.26
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-4.46
eij[("ASP","ASP")]=0.84
eij[("ASP","GLU")]=eij[("GLU","ASP")]=1.97
eij[("ASP","PHE")]=eij[("PHE","ASP")]=-0.92
eij[("ASP","GLY")]=eij[("GLY","ASP")]=0.88
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-1.07
eij[("ASP","ILE")]=eij[("ILE","ASP")]=0.68
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-1.93
eij[("ASP","LEU")]=eij[("LEU","ASP")]=0.23
eij[("ASP","MET")]=eij[("MET","ASP")]=0.61
eij[("ASP","ASN")]=eij[("ASN","ASP")]=0.32
eij[("ASP","PRO")]=eij[("PRO","ASP")]=3.31
eij[("ASP","GLN")]=eij[("GLN","ASP")]=2.67
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-2.02
eij[("ASP","SER")]=eij[("SER","ASP")]=0.91
eij[("ASP","THR")]=eij[("THR","ASP")]=-0.65
eij[("ASP","VAL")]=eij[("VAL","ASP")]=0.94
eij[("ASP","TRP")]=eij[("TRP","ASP")]=-0.71
eij[("ASP","TYR")]=eij[("TYR","ASP")]=0.9
eij[("GLU","GLU")]=1.45
eij[("GLU","PHE")]=eij[("PHE","GLU")]=0.94
eij[("GLU","GLY")]=eij[("GLY","GLU")]=1.31
eij[("GLU","HIS")]=eij[("HIS","GLU")]=0.61
eij[("GLU","ILE")]=eij[("ILE","GLU")]=1.3
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-2.51
eij[("GLU","LEU")]=eij[("LEU","GLU")]=1.14
eij[("GLU","MET")]=eij[("MET","GLU")]=2.53
eij[("GLU","ASN")]=eij[("ASN","GLU")]=0.2
eij[("GLU","PRO")]=eij[("PRO","GLU")]=1.44
eij[("GLU","GLN")]=eij[("GLN","GLU")]=0.1
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-3.13
eij[("GLU","SER")]=eij[("SER","GLU")]=0.81
eij[("GLU","THR")]=eij[("THR","GLU")]=1.54
eij[("GLU","VAL")]=eij[("VAL","GLU")]=0.12
eij[("GLU","TRP")]=eij[("TRP","GLU")]=-1.07
eij[("GLU","TYR")]=eij[("TYR","GLU")]=1.29
eij[("PHE","PHE")]=-11.25
eij[("PHE","GLY")]=eij[("GLY","PHE")]=0.35
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-3.57
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-5.88
eij[("PHE","LYS")]=eij[("LYS","PHE")]=-0.82
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-8.59
eij[("PHE","MET")]=eij[("MET","PHE")]=-5.34
eij[("PHE","ASN")]=eij[("ASN","PHE")]=0.73
eij[("PHE","PRO")]=eij[("PRO","PHE")]=0.32
eij[("PHE","GLN")]=eij[("GLN","PHE")]=0.77
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-0.4
eij[("PHE","SER")]=eij[("SER","PHE")]=-2.22
eij[("PHE","THR")]=eij[("THR","PHE")]=0.11
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-7.05
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-7.09
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-8.8
eij[("GLY","GLY")]=-0.2
eij[("GLY","HIS")]=eij[("HIS","GLY")]=1.09
eij[("GLY","ILE")]=eij[("ILE","GLY")]=-0.65
eij[("GLY","LYS")]=eij[("LYS","GLY")]=-0.16
eij[("GLY","LEU")]=eij[("LEU","GLY")]=-0.55
eij[("GLY","MET")]=eij[("MET","GLY")]=-0.52
eij[("GLY","ASN")]=eij[("ASN","GLY")]=-0.32
eij[("GLY","PRO")]=eij[("PRO","GLY")]=2.25
eij[("GLY","GLN")]=eij[("GLN","GLY")]=1.11
eij[("GLY","ARG")]=eij[("ARG","GLY")]=0.84
eij[("GLY","SER")]=eij[("SER","GLY")]=0.71
eij[("GLY","THR")]=eij[("THR","GLY")]=0.59
eij[("GLY","VAL")]=eij[("VAL","GLY")]=-0.38
eij[("GLY","TRP")]=eij[("TRP","GLY")]=1.69
eij[("GLY","TYR")]=eij[("TYR","GLY")]=-1.9
eij[("HIS","HIS")]=1.97
eij[("HIS","ILE")]=eij[("ILE","HIS")]=-0.71
eij[("HIS","LYS")]=eij[("LYS","HIS")]=2.89
eij[("HIS","LEU")]=eij[("LEU","HIS")]=-0.86
eij[("HIS","MET")]=eij[("MET","HIS")]=-0.75
eij[("HIS","ASN")]=eij[("ASN","HIS")]=1.84
eij[("HIS","PRO")]=eij[("PRO","HIS")]=0.35
eij[("HIS","GLN")]=eij[("GLN","HIS")]=2.64
eij[("HIS","ARG")]=eij[("ARG","HIS")]=2.05
eij[("HIS","SER")]=eij[("SER","HIS")]=0.82
eij[("HIS","THR")]=eij[("THR","HIS")]=-0.01
eij[("HIS","VAL")]=eij[("VAL","HIS")]=0.27
eij[("HIS","TRP")]=eij[("TRP","HIS")]=-7.58
eij[("HIS","TYR")]=eij[("TYR","HIS")]=-3.2
eij[("ILE","ILE")]=-6.74
eij[("ILE","LYS")]=eij[("LYS","ILE")]=-0.01
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-9.01
eij[("ILE","MET")]=eij[("MET","ILE")]=-3.62
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-0.07
eij[("ILE","PRO")]=eij[("PRO","ILE")]=0.12
eij[("ILE","GLN")]=eij[("GLN","ILE")]=-0.18
eij[("ILE","ARG")]=eij[("ARG","ILE")]=0.19
eij[("ILE","SER")]=eij[("SER","ILE")]=-0.15
eij[("ILE","THR")]=eij[("THR","ILE")]=0.63
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-6.54
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-3.78
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-5.26
eij[("LYS","LYS")]=1.24
eij[("LYS","LEU")]=eij[("LEU","LYS")]=0.49
eij[("LYS","MET")]=eij[("MET","LYS")]=1.61
eij[("LYS","ASN")]=eij[("ASN","LYS")]=1.12
eij[("LYS","PRO")]=eij[("PRO","LYS")]=0.51
eij[("LYS","GLN")]=eij[("GLN","LYS")]=0.43
eij[("LYS","ARG")]=eij[("ARG","LYS")]=2.34
eij[("LYS","SER")]=eij[("SER","LYS")]=0.19
eij[("LYS","THR")]=eij[("THR","LYS")]=-1.11
eij[("LYS","VAL")]=eij[("VAL","LYS")]=0.19
eij[("LYS","TRP")]=eij[("TRP","LYS")]=0.02
eij[("LYS","TYR")]=eij[("TYR","LYS")]=-1.19
eij[("LEU","LEU")]=-6.37
eij[("LEU","MET")]=eij[("MET","LEU")]=-2.88
eij[("LEU","ASN")]=eij[("ASN","LEU")]=0.97
eij[("LEU","PRO")]=eij[("PRO","LEU")]=1.81
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-0.58
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-0.6
eij[("LEU","SER")]=eij[("SER","LEU")]=-0.41
eij[("LEU","THR")]=eij[("THR","LEU")]=0.72
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-5.43
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-8.31
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-4.9
eij[("MET","MET")]=-6.49
eij[("MET","ASN")]=eij[("ASN","MET")]=0.21
eij[("MET","PRO")]=eij[("PRO","MET")]=0.75
eij[("MET","GLN")]=eij[("GLN","MET")]=1.9
eij[("MET","ARG")]=eij[("ARG","MET")]=2.09
eij[("MET","SER")]=eij[("SER","MET")]=1.39
eij[("MET","THR")]=eij[("THR","MET")]=0.63
eij[("MET","VAL")]=eij[("VAL","MET")]=-2.59
eij[("MET","TRP")]=eij[("TRP","MET")]=-6.88
eij[("MET","TYR")]=eij[("TYR","MET")]=-9.73
eij[("ASN","ASN")]=0.61
eij[("ASN","PRO")]=eij[("PRO","ASN")]=1.15
eij[("ASN","GLN")]=eij[("GLN","ASN")]=1.28
eij[("ASN","ARG")]=eij[("ARG","ASN")]=1.08
eij[("ASN","SER")]=eij[("SER","ASN")]=0.29
eij[("ASN","THR")]=eij[("THR","ASN")]=0.46
eij[("ASN","VAL")]=eij[("VAL","ASN")]=0.93
eij[("ASN","TRP")]=eij[("TRP","ASN")]=-0.74
eij[("ASN","TYR")]=eij[("TYR","ASN")]=0.93
eij[("PRO","PRO")]=-0.42
eij[("PRO","GLN")]=eij[("GLN","PRO")]=2.97
eij[("PRO","ARG")]=eij[("ARG","PRO")]=1.06
eij[("PRO","SER")]=eij[("SER","PRO")]=1.12
eij[("PRO","THR")]=eij[("THR","PRO")]=1.65
eij[("PRO","VAL")]=eij[("VAL","PRO")]=0.38
eij[("PRO","TRP")]=eij[("TRP","PRO")]=-2.06
eij[("PRO","TYR")]=eij[("TYR","PRO")]=-2.09
eij[("GLN","GLN")]=-1.54
eij[("GLN","ARG")]=eij[("ARG","GLN")]=0.91
eij[("GLN","SER")]=eij[("SER","GLN")]=0.85
eij[("GLN","THR")]=eij[("THR","GLN")]=-0.07
eij[("GLN","VAL")]=eij[("VAL","GLN")]=-1.91
eij[("GLN","TRP")]=eij[("TRP","GLN")]=-0.76
eij[("GLN","TYR")]=eij[("TYR","GLN")]=0.01
eij[("ARG","ARG")]=0.21
eij[("ARG","SER")]=eij[("SER","ARG")]=0.95
eij[("ARG","THR")]=eij[("THR","ARG")]=0.98
eij[("ARG","VAL")]=eij[("VAL","ARG")]=0.08
eij[("ARG","TRP")]=eij[("TRP","ARG")]=-5.89
eij[("ARG","TYR")]=eij[("TYR","ARG")]=0.36
eij[("SER","SER")]=-0.48
eij[("SER","THR")]=eij[("THR","SER")]=-0.06
eij[("SER","VAL")]=eij[("VAL","SER")]=0.13
eij[("SER","TRP")]=eij[("TRP","SER")]=-3.03
eij[("SER","TYR")]=eij[("TYR","SER")]=-0.82
eij[("THR","THR")]=-0.96
eij[("THR","VAL")]=eij[("VAL","THR")]=1.14
eij[("THR","TRP")]=eij[("TRP","THR")]=-0.65
eij[("THR","TYR")]=eij[("TYR","THR")]=-0.37
eij[("VAL","VAL")]=-4.82
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-2.13
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-3.59
eij[("TRP","TRP")]=-1.73
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-12.39
eij[("TYR","TYR")]=-2.68

def scoreEC_BB_NB(pdbid, ctype):

    if (ctype == 'CA'):
        file1 = open('./raw_graphv2/'+pdbid+'_CA_output_bb_nb.txt', 'r')
    elif (ctype == 'CB'):
        file1 = open('./raw_graphv2/'+pdbid+'_CB_output_bb.txt', 'r')
    
    Lines = file1.readlines()
    total_lines=len(Lines)

    if (ctype == 'CA'):
        with open('./raw_graphv2/'+pdbid+'_CA_output_bb_nb.txt', 'w') as f:
            print("")
    elif (ctype == 'CB'):
        with open('./raw_graphv2/'+pdbid+'_CB_output_bb.txt', 'w') as f:
            print("")

    Ec = 0
    Etot=0
    count = 0
    resPid = Lines[0].split(',')[2]
    num_lines=0
    last_line=0

    # Strips the newline character
    for line in Lines:
        resA=line.split(',')[1]
        resAid=line.split(',')[2]
        resB=line.split(',')[4]
        dist=line.split(',')[6]
        num_lines=num_lines+1

        Eij=-1*(eij[(resA,resB)])/100
    
        line=line.strip("\n")+","+str(round(Eij,2))
        if (ctype == 'CA'):
            with open('./raw_graphv2/'+pdbid+'_CA_output_bb_nb.txt', 'a') as f:
                print(line, file=f)
        if (ctype == 'CB'):
            with open('./raw_graphv2/'+pdbid+'_CB_output_bb.txt', 'a') as f:
                print(line, file=f)


        if (num_lines==total_lines):
            last_line=1
            Ec=Ec+Eij
        else:
            last_line=0

        if (resAid != resPid or last_line==1):
            Etot=Etot+Ec
            #print("resPid:"+resPid+" Etot:"+str(round(Etot/100,2))+" Ec:"+str(round(Ec,2))+" lines:"+str(num_lines)+"/"+str(total_lines))
            Ec=0
            resPid=resAid
    
    
        Ec=Ec+Eij
        
    print("scoreEC_BB_NB:",str(round(Etot,2))," ctype:",ctype)
    return round(Etot,2)
