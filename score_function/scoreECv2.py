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
eij[("ALA","ALA")]=-0.2
eij[("ALA","CYS")]=eij[("CYS","ALA")]=-0.44
eij[("ALA","ASP")]=eij[("ASP","ALA")]=0.16
eij[("ALA","GLU")]=eij[("GLU","ALA")]=0.26
eij[("ALA","PHE")]=eij[("PHE","ALA")]=-0.46
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-0.26
eij[("ALA","HIS")]=eij[("HIS","ALA")]=0.5
eij[("ALA","ILE")]=eij[("ILE","ALA")]=-0.57
eij[("ALA","LYS")]=eij[("LYS","ALA")]=0.1
eij[("ALA","LEU")]=eij[("LEU","ALA")]=-0.36
eij[("ALA","MET")]=eij[("MET","ALA")]=-0.22
eij[("ALA","ASN")]=eij[("ASN","ALA")]=0.07
eij[("ALA","PRO")]=eij[("PRO","ALA")]=0.14
eij[("ALA","GLN")]=eij[("GLN","ALA")]=0.01
eij[("ALA","ARG")]=eij[("ARG","ALA")]=0.2
eij[("ALA","SER")]=eij[("SER","ALA")]=-0.09
eij[("ALA","THR")]=eij[("THR","ALA")]=-0.05
eij[("ALA","VAL")]=eij[("VAL","ALA")]=-0.42
eij[("ALA","TRP")]=eij[("TRP","ALA")]=0.05
eij[("ALA","TYR")]=eij[("TYR","ALA")]=-0.5
eij[("CYS","CYS")]=-2.99
eij[("CYS","ASP")]=eij[("ASP","CYS")]=0.21
eij[("CYS","GLU")]=eij[("GLU","CYS")]=0.19
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-0.88
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-0.34
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-1.11
eij[("CYS","ILE")]=eij[("ILE","CYS")]=-0.36
eij[("CYS","LYS")]=eij[("LYS","CYS")]=-0.09
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-0.53
eij[("CYS","MET")]=eij[("MET","CYS")]=-0.43
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-0.52
eij[("CYS","PRO")]=eij[("PRO","CYS")]=-0.14
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-0.43
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-0.24
eij[("CYS","SER")]=eij[("SER","CYS")]=0.13
eij[("CYS","THR")]=eij[("THR","CYS")]=-0.22
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-0.62
eij[("CYS","TRP")]=eij[("TRP","CYS")]=0.24
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-0.79
eij[("ASP","ASP")]=0.17
eij[("ASP","GLU")]=eij[("GLU","ASP")]=0.55
eij[("ASP","PHE")]=eij[("PHE","ASP")]=0.38
eij[("ASP","GLY")]=eij[("GLY","ASP")]=0.35
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-0.23
eij[("ASP","ILE")]=eij[("ILE","ASP")]=0.44
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-0.39
eij[("ASP","LEU")]=eij[("LEU","ASP")]=0.28
eij[("ASP","MET")]=eij[("MET","ASP")]=0.35
eij[("ASP","ASN")]=eij[("ASN","ASP")]=-0.02
eij[("ASP","PRO")]=eij[("PRO","ASP")]=1.03
eij[("ASP","GLN")]=eij[("GLN","ASP")]=0.49
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-0.37
eij[("ASP","SER")]=eij[("SER","ASP")]=0.19
eij[("ASP","THR")]=eij[("THR","ASP")]=-0.12
eij[("ASP","VAL")]=eij[("VAL","ASP")]=0.69
eij[("ASP","TRP")]=eij[("TRP","ASP")]=0.04
eij[("ASP","TYR")]=eij[("TYR","ASP")]=0.43
eij[("GLU","GLU")]=0.6
eij[("GLU","PHE")]=eij[("PHE","GLU")]=0.55
eij[("GLU","GLY")]=eij[("GLY","GLU")]=0.65
eij[("GLU","HIS")]=eij[("HIS","GLU")]=0.18
eij[("GLU","ILE")]=eij[("ILE","GLU")]=0.37
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-0.47
eij[("GLU","LEU")]=eij[("LEU","GLU")]=0.33
eij[("GLU","MET")]=eij[("MET","GLU")]=0.29
eij[("GLU","ASN")]=eij[("ASN","GLU")]=0.01
eij[("GLU","PRO")]=eij[("PRO","GLU")]=0.69
eij[("GLU","GLN")]=eij[("GLN","GLU")]=0.04
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-0.52
eij[("GLU","SER")]=eij[("SER","GLU")]=0.18
eij[("GLU","THR")]=eij[("THR","GLU")]=0.37
eij[("GLU","VAL")]=eij[("VAL","GLU")]=0.39
eij[("GLU","TRP")]=eij[("TRP","GLU")]=0.03
eij[("GLU","TYR")]=eij[("TYR","GLU")]=0.17
eij[("PHE","PHE")]=-0.94
eij[("PHE","GLY")]=eij[("GLY","PHE")]=0.17
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-0.4
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-0.88
eij[("PHE","LYS")]=eij[("LYS","PHE")]=0.01
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-1.08
eij[("PHE","MET")]=eij[("MET","PHE")]=-0.78
eij[("PHE","ASN")]=eij[("ASN","PHE")]=0.22
eij[("PHE","PRO")]=eij[("PRO","PHE")]=0.2
eij[("PHE","GLN")]=eij[("GLN","PHE")]=0.26
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-0.19
eij[("PHE","SER")]=eij[("SER","PHE")]=-0.22
eij[("PHE","THR")]=eij[("THR","PHE")]=0.02
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-1.15
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-0.6
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-0.88
eij[("GLY","GLY")]=-0.12
eij[("GLY","HIS")]=eij[("HIS","GLY")]=0.18
eij[("GLY","ILE")]=eij[("ILE","GLY")]=0.24
eij[("GLY","LYS")]=eij[("LYS","GLY")]=0.19
eij[("GLY","LEU")]=eij[("LEU","GLY")]=0.24
eij[("GLY","MET")]=eij[("MET","GLY")]=0.02
eij[("GLY","ASN")]=eij[("ASN","GLY")]=-0.04
eij[("GLY","PRO")]=eij[("PRO","GLY")]=0.6
eij[("GLY","GLN")]=eij[("GLN","GLY")]=0.46
eij[("GLY","ARG")]=eij[("ARG","GLY")]=0.5
eij[("GLY","SER")]=eij[("SER","GLY")]=0.28
eij[("GLY","THR")]=eij[("THR","GLY")]=0.28
eij[("GLY","VAL")]=eij[("VAL","GLY")]=0.27
eij[("GLY","TRP")]=eij[("TRP","GLY")]=0.51
eij[("GLY","TYR")]=eij[("TYR","GLY")]=-0.35
eij[("HIS","HIS")]=0.42
eij[("HIS","ILE")]=eij[("ILE","HIS")]=0
eij[("HIS","LYS")]=eij[("LYS","HIS")]=0.79
eij[("HIS","LEU")]=eij[("LEU","HIS")]=-0.24
eij[("HIS","MET")]=eij[("MET","HIS")]=-0.07
eij[("HIS","ASN")]=eij[("ASN","HIS")]=0.2
eij[("HIS","PRO")]=eij[("PRO","HIS")]=0.25
eij[("HIS","GLN")]=eij[("GLN","HIS")]=0.69
eij[("HIS","ARG")]=eij[("ARG","HIS")]=0.24
eij[("HIS","SER")]=eij[("SER","HIS")]=0.21
eij[("HIS","THR")]=eij[("THR","HIS")]=0.11
eij[("HIS","VAL")]=eij[("VAL","HIS")]=0.16
eij[("HIS","TRP")]=eij[("TRP","HIS")]=-0.85
eij[("HIS","TYR")]=eij[("TYR","HIS")]=-0.26
eij[("ILE","ILE")]=-1.16
eij[("ILE","LYS")]=eij[("LYS","ILE")]=0.15
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-1.25
eij[("ILE","MET")]=eij[("MET","ILE")]=-0.58
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-0.09
eij[("ILE","PRO")]=eij[("PRO","ILE")]=0.36
eij[("ILE","GLN")]=eij[("GLN","ILE")]=-0.08
eij[("ILE","ARG")]=eij[("ARG","ILE")]=0.14
eij[("ILE","SER")]=eij[("SER","ILE")]=0.32
eij[("ILE","THR")]=eij[("THR","ILE")]=-0.27
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-1.06
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-0.68
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-0.85
eij[("LYS","LYS")]=0.42
eij[("LYS","LEU")]=eij[("LEU","LYS")]=0.13
eij[("LYS","MET")]=eij[("MET","LYS")]=0.48
eij[("LYS","ASN")]=eij[("ASN","LYS")]=0.26
eij[("LYS","PRO")]=eij[("PRO","LYS")]=0.5
eij[("LYS","GLN")]=eij[("GLN","LYS")]=0.15
eij[("LYS","ARG")]=eij[("ARG","LYS")]=0.53
eij[("LYS","SER")]=eij[("SER","LYS")]=0.1
eij[("LYS","THR")]=eij[("THR","LYS")]=-0.19
eij[("LYS","VAL")]=eij[("VAL","LYS")]=0.1
eij[("LYS","TRP")]=eij[("TRP","LYS")]=0.1
eij[("LYS","TYR")]=eij[("TYR","LYS")]=0.04
eij[("LEU","LEU")]=-1.10
eij[("LEU","MET")]=eij[("MET","LEU")]=-0.5
eij[("LEU","ASN")]=eij[("ASN","LEU")]=0.21
eij[("LEU","PRO")]=eij[("PRO","LEU")]=0.42
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-0.01
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-0.07
eij[("LEU","SER")]=eij[("SER","LEU")]=0.17
eij[("LEU","THR")]=eij[("THR","LEU")]=0.07
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-0.97
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-0.95
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-0.63
eij[("MET","MET")]=-0.74
eij[("MET","ASN")]=eij[("ASN","MET")]=0.32
eij[("MET","PRO")]=eij[("PRO","MET")]=0.01
eij[("MET","GLN")]=eij[("GLN","MET")]=0.26
eij[("MET","ARG")]=eij[("ARG","MET")]=0.15
eij[("MET","SER")]=eij[("SER","MET")]=0.48
eij[("MET","THR")]=eij[("THR","MET")]=0.16
eij[("MET","VAL")]=eij[("VAL","MET")]=-0.73
eij[("MET","TRP")]=eij[("TRP","MET")]=-0.56
eij[("MET","TYR")]=eij[("TYR","MET")]=-1.02
eij[("ASN","ASN")]=0.14
eij[("ASN","PRO")]=eij[("PRO","ASN")]=0.27
eij[("ASN","GLN")]=eij[("GLN","ASN")]=0.37
eij[("ASN","ARG")]=eij[("ARG","ASN")]=0.13
eij[("ASN","SER")]=eij[("SER","ASN")]=0.15
eij[("ASN","THR")]=eij[("THR","ASN")]=0.1
eij[("ASN","VAL")]=eij[("VAL","ASN")]=0.4
eij[("ASN","TRP")]=eij[("TRP","ASN")]=-0.12
eij[("ASN","TYR")]=eij[("TYR","ASN")]=0.32
eij[("PRO","PRO")]=0.27
eij[("PRO","GLN")]=eij[("GLN","PRO")]=1.02
eij[("PRO","ARG")]=eij[("ARG","PRO")]=0.47
eij[("PRO","SER")]=eij[("SER","PRO")]=0.54
eij[("PRO","THR")]=eij[("THR","PRO")]=0.88
eij[("PRO","VAL")]=eij[("VAL","PRO")]=-0.02
eij[("PRO","TRP")]=eij[("TRP","PRO")]=-0.37
eij[("PRO","TYR")]=eij[("TYR","PRO")]=-0.12
eij[("GLN","GLN")]=-0.12
eij[("GLN","ARG")]=eij[("ARG","GLN")]=0.24
eij[("GLN","SER")]=eij[("SER","GLN")]=0.29
eij[("GLN","THR")]=eij[("THR","GLN")]=0.04
eij[("GLN","VAL")]=eij[("VAL","GLN")]=-0.11
eij[("GLN","TRP")]=eij[("TRP","GLN")]=0.18
eij[("GLN","TYR")]=eij[("TYR","GLN")]=0.11
eij[("ARG","ARG")]=0.17
eij[("ARG","SER")]=eij[("SER","ARG")]=0.27
eij[("ARG","THR")]=eij[("THR","ARG")]=0.45
eij[("ARG","VAL")]=eij[("VAL","ARG")]=0.01
eij[("ARG","TRP")]=eij[("TRP","ARG")]=-0.73
eij[("ARG","TYR")]=eij[("TYR","ARG")]=0.01
eij[("SER","SER")]=-0.06
eij[("SER","THR")]=eij[("THR","SER")]=0.08
eij[("SER","VAL")]=eij[("VAL","SER")]=0.12
eij[("SER","TRP")]=eij[("TRP","SER")]=-0.22
eij[("SER","TYR")]=eij[("TYR","SER")]=-0.14
eij[("THR","THR")]=-0.03
eij[("THR","VAL")]=eij[("VAL","THR")]=-0.01
eij[("THR","TRP")]=eij[("TRP","THR")]=0.11
eij[("THR","TYR")]=eij[("TYR","THR")]=-0.32
eij[("VAL","VAL")]=-0.89
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-0.56
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-0.71
eij[("TRP","TRP")]=-0.05
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-1.41
eij[("TYR","TYR")]=-0.76


def scoreEC(pdbid, ctype):


    if (ctype == 'CA'):
        file1 = open('./raw_graphv2/'+pdbid+'_CA_output.txt', 'r')
    elif (ctype == 'CB'):
        file1 = open('./raw_graphv2/'+pdbid+'_CB_output.txt', 'r')

    Lines = file1.readlines()
    total_lines=len(Lines)

    if (ctype == 'CA'):
        with open('./raw_graphv2/'+pdbid+'_CA_output.txt', 'w') as f:
            print("")
    elif (ctype == 'CB'):
        with open('./raw_graphv2/'+pdbid+'_CB_output.txt', 'w') as f:
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

        Eij=-1*(eij[(resA,resB)])
    
        line=line.strip("\n")+","+str(round(Eij,2))
        if (ctype == 'CA'):
            with open('./raw_graphv2/'+pdbid+'_CA_output.txt', 'a') as f:
                print(line, file=f)
        if (ctype == 'CB'):
            with open('./raw_graphv2/'+pdbid+'_CB_output.txt', 'a') as f:
                print(line, file=f)

        if (num_lines==total_lines):
            last_line=1
            Ec=Ec+Eij
        else:
            last_line=0   
	
        if (resAid != resPid or last_line==1):
            Etot=Etot+Ec
            #print("resPid:"+resPid+" Etot:"+str(round(Etot,2))+" Ec:"+str(round(Ec,2))+" lines:"+str(num_lines)+"/"+str(total_lines))
            Ec=0
            resPid=resAid
    
        Ec=Ec+Eij
        
    print("scoreEC:",str(round(Etot,2))," ctype:",ctype)
    return round(Etot,2)
