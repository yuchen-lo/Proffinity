import argparse
import sys
from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt, log

#################
#Miyazawa et al. (1996)
#Eij=eij+err-eir-ejr
eij={}
eij[("CYS","CYS")]=-1.79
eij[("CYS","MET")]=eij[("MET","CYS")]=-1.23
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-0.98
eij[("CYS","ILE")]=eij[("ILE","CYS")]=-0.48
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-0.69
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-0.94
eij[("CYS","TRP")]=eij[("TRP","CYS")]=-0.3
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-0.96
eij[("CYS","ALA")]=eij[("ALA","CYS")]=-0.3
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-0.42
eij[("CYS","THR")]=eij[("THR","CYS")]=-0.38
eij[("CYS","SER")]=eij[("SER","CYS")]=-0.2
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-0.49
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-0.32
eij[("CYS","ASP")]=eij[("ASP","CYS")]=0.04
eij[("CYS","GLU")]=eij[("GLU","CYS")]=0.55
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-0.82
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-0.4
eij[("CYS","LYS")]=eij[("LYS","CYS")]=0
eij[("CYS","PRO")]=eij[("PRO","CYS")]=0.07
eij[("MET","MET")]=0.36
eij[("MET","PHE")]=eij[("PHE","MET")]=-1.03
eij[("MET","ILE")]=eij[("ILE","MET")]=-0.41
eij[("MET","LEU")]=eij[("LEU","MET")]=-0.31
eij[("MET","VAL")]=eij[("VAL","MET")]=-0.94
eij[("MET","TRP")]=eij[("TRP","MET")]=-0.07
eij[("MET","TYR")]=eij[("TYR","MET")]=-1.10
eij[("MET","ALA")]=eij[("ALA","MET")]=0.05
eij[("MET","GLY")]=eij[("GLY","MET")]=0
eij[("MET","THR")]=eij[("THR","MET")]=0.06
eij[("MET","SER")]=eij[("SER","MET")]=-0.47
eij[("MET","ASN")]=eij[("ASN","MET")]=-0.54
eij[("MET","GLN")]=eij[("GLN","MET")]=0.31
eij[("MET","ASP")]=eij[("ASP","MET")]=0.02
eij[("MET","GLU")]=eij[("GLU","MET")]=1.07
eij[("MET","HIS")]=eij[("HIS","MET")]=-0.35
eij[("MET","ARG")]=eij[("ARG","MET")]=-0.43
eij[("MET","LYS")]=eij[("LYS","MET")]=0.55
eij[("MET","PRO")]=eij[("PRO","MET")]=-0.25
eij[("PHE","PHE")]=-0.61
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-0.66
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-1.02
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-0.78
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-0.89
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-0.82
eij[("PHE","ALA")]=eij[("ALA","PHE")]=-0.05
eij[("PHE","GLY")]=eij[("GLY","PHE")]=0.21
eij[("PHE","THR")]=eij[("THR","PHE")]=-0.19
eij[("PHE","SER")]=eij[("SER","PHE")]=0.14
eij[("PHE","ASN")]=eij[("ASN","PHE")]=0.1
eij[("PHE","GLN")]=eij[("GLN","PHE")]=-0.02
eij[("PHE","ASP")]=eij[("ASP","PHE")]=0.19
eij[("PHE","GLU")]=eij[("GLU","PHE")]=0.2
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-0.75
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-0.22
eij[("PHE","LYS")]=eij[("LYS","PHE")]=-0.17
eij[("PHE","PRO")]=eij[("PRO","PHE")]=-0.43
eij[("ILE","ILE")]=-0.71
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-1.04
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-0.98
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-0.89
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-0.87
eij[("ILE","ALA")]=eij[("ALA","ILE")]=-0.64
eij[("ILE","GLY")]=eij[("GLY","ILE")]=0.4
eij[("ILE","THR")]=eij[("THR","ILE")]=-0.29
eij[("ILE","SER")]=eij[("SER","ILE")]=-0.13
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-0.39
eij[("ILE","GLN")]=eij[("GLN","ILE")]=0.39
eij[("ILE","ASP")]=eij[("ASP","ILE")]=-0.2
eij[("ILE","GLU")]=eij[("GLU","ILE")]=0.04
eij[("ILE","HIS")]=eij[("HIS","ILE")]=-0.52
eij[("ILE","ARG")]=eij[("ARG","ILE")]=-0.08
eij[("ILE","LYS")]=eij[("LYS","ILE")]=-0.26
eij[("ILE","PRO")]=eij[("PRO","ILE")]=0.25
eij[("LEU","LEU")]=-1.14
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-1.03
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-0.97
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-0.6
eij[("LEU","ALA")]=eij[("ALA","LEU")]=-0.57
eij[("LEU","GLY")]=eij[("GLY","LEU")]=-0.08
eij[("LEU","THR")]=eij[("THR","LEU")]=-0.39
eij[("LEU","SER")]=eij[("SER","LEU")]=-0.07
eij[("LEU","ASN")]=eij[("ASN","LEU")]=-0.13
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-0.1
eij[("LEU","ASP")]=eij[("ASP","LEU")]=-0.05
eij[("LEU","GLU")]=eij[("GLU","LEU")]=0.5
eij[("LEU","HIS")]=eij[("HIS","LEU")]=-0.36
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-0.1
eij[("LEU","LYS")]=eij[("LYS","LEU")]=0.1
eij[("LEU","PRO")]=eij[("PRO","LEU")]=0.09
eij[("VAL","VAL")]=-1.15
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-0.6
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-0.7
eij[("VAL","ALA")]=eij[("ALA","VAL")]=-0.6
eij[("VAL","GLY")]=eij[("GLY","VAL")]=-0.2
eij[("VAL","THR")]=eij[("THR","VAL")]=0.06
eij[("VAL","SER")]=eij[("SER","VAL")]=-0.31
eij[("VAL","ASN")]=eij[("ASN","VAL")]=-0.09
eij[("VAL","GLN")]=eij[("GLN","VAL")]=-0.24
eij[("VAL","ASP")]=eij[("ASP","VAL")]=-0.02
eij[("VAL","GLU")]=eij[("GLU","VAL")]=0.25
eij[("VAL","HIS")]=eij[("HIS","VAL")]=-0.35
eij[("VAL","ARG")]=eij[("ARG","VAL")]=-0.48
eij[("VAL","LYS")]=eij[("LYS","VAL")]=-0.08
eij[("VAL","PRO")]=eij[("PRO","VAL")]=-0.08
eij[("TRP","TRP")]=0.02
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-0.99
eij[("TRP","ALA")]=eij[("ALA","TRP")]=-0.08
eij[("TRP","GLY")]=eij[("GLY","TRP")]=-0.14
eij[("TRP","THR")]=eij[("THR","TRP")]=0.07
eij[("TRP","SER")]=eij[("SER","TRP")]=-0.2
eij[("TRP","ASN")]=eij[("ASN","TRP")]=0.4
eij[("TRP","GLN")]=eij[("GLN","TRP")]=-0.68
eij[("TRP","ASP")]=eij[("ASP","TRP")]=0.32
eij[("TRP","GLU")]=eij[("GLU","TRP")]=0.24
eij[("TRP","HIS")]=eij[("HIS","TRP")]=-0.41
eij[("TRP","ARG")]=eij[("ARG","TRP")]=-0.78
eij[("TRP","LYS")]=eij[("LYS","TRP")]=-0.3
eij[("TRP","PRO")]=eij[("PRO","TRP")]=-0.44
eij[("TYR","TYR")]=0.35
eij[("TYR","ALA")]=eij[("ALA","TYR")]=-0.37
eij[("TYR","GLY")]=eij[("GLY","TYR")]=-0.32
eij[("TYR","THR")]=eij[("THR","TYR")]=-0.23
eij[("TYR","SER")]=eij[("SER","TYR")]=0.25
eij[("TYR","ASN")]=eij[("ASN","TYR")]=-0.39
eij[("TYR","GLN")]=eij[("GLN","TYR")]=-0.74
eij[("TYR","ASP")]=eij[("ASP","TYR")]=0.22
eij[("TYR","GLU")]=eij[("GLU","TYR")]=0.11
eij[("TYR","HIS")]=eij[("HIS","TYR")]=-0.67
eij[("TYR","ARG")]=eij[("ARG","TYR")]=0.21
eij[("TYR","LYS")]=eij[("LYS","TYR")]=-0.2
eij[("TYR","PRO")]=eij[("PRO","TYR")]=-0.45
eij[("ALA","ALA")]=-0.08
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-0.09
eij[("ALA","THR")]=eij[("THR","ALA")]=-0.22
eij[("ALA","SER")]=eij[("SER","ALA")]=-0.01
eij[("ALA","ASN")]=eij[("ASN","ALA")]=-0.11
eij[("ALA","GLN")]=eij[("GLN","ALA")]=-0.14
eij[("ALA","ASP")]=eij[("ASP","ALA")]=0.03
eij[("ALA","GLU")]=eij[("GLU","ALA")]=0.1
eij[("ALA","HIS")]=eij[("HIS","ALA")]=-0.15
eij[("ALA","ARG")]=eij[("ARG","ALA")]=0.07
eij[("ALA","LYS")]=eij[("LYS","ALA")]=0
eij[("ALA","PRO")]=eij[("PRO","ALA")]=0.41
eij[("GLY","GLY")]=0.04
eij[("GLY","THR")]=eij[("THR","GLY")]=0.13
eij[("GLY","SER")]=eij[("SER","GLY")]=-0.04
eij[("GLY","ASN")]=eij[("ASN","GLY")]=0.12
eij[("GLY","GLN")]=eij[("GLN","GLY")]=-0.18
eij[("GLY","ASP")]=eij[("ASP","GLY")]=0.4
eij[("GLY","GLU")]=eij[("GLU","GLY")]=-0.06
eij[("GLY","HIS")]=eij[("HIS","GLY")]=0
eij[("GLY","ARG")]=eij[("ARG","GLY")]=-0.15
eij[("GLY","LYS")]=eij[("LYS","GLY")]=0.1
eij[("GLY","PRO")]=eij[("PRO","GLY")]=0.4
eij[("THR","THR")]=0.26
eij[("THR","SER")]=eij[("SER","THR")]=0.05
eij[("THR","ASN")]=eij[("ASN","THR")]=-0.17
eij[("THR","GLN")]=eij[("GLN","THR")]=-0.27
eij[("THR","ASP")]=eij[("ASP","THR")]=0.15
eij[("THR","GLU")]=eij[("GLU","THR")]=-0.03
eij[("THR","HIS")]=eij[("HIS","THR")]=-0.27
eij[("THR","ARG")]=eij[("ARG","THR")]=-0.17
eij[("THR","LYS")]=eij[("LYS","THR")]=0.09
eij[("THR","PRO")]=eij[("PRO","THR")]=0.36
eij[("SER","SER")]=-0.13
eij[("SER","ASN")]=eij[("ASN","SER")]=0.4
eij[("SER","GLN")]=eij[("GLN","SER")]=0.37
eij[("SER","ASP")]=eij[("ASP","SER")]=0.3
eij[("SER","GLU")]=eij[("GLU","SER")]=-0.09
eij[("SER","HIS")]=eij[("HIS","SER")]=-0.59
eij[("SER","ARG")]=eij[("ARG","SER")]=0.61
eij[("SER","LYS")]=eij[("LYS","SER")]=0.18
eij[("SER","PRO")]=eij[("PRO","SER")]=0.44
eij[("ASN","ASN")]=-0.08
eij[("ASN","GLN")]=eij[("GLN","ASN")]=-0.05
eij[("ASN","ASP")]=eij[("ASP","ASN")]=0.62
eij[("ASN","GLU")]=eij[("GLU","ASN")]=0.46
eij[("ASN","HIS")]=eij[("HIS","ASN")]=0.05
eij[("ASN","ARG")]=eij[("ARG","ASN")]=0.62
eij[("ASN","LYS")]=eij[("LYS","ASN")]=0.04
eij[("ASN","PRO")]=eij[("PRO","ASN")]=-0.21
eij[("GLN","GLN")]=-0.86
eij[("GLN","ASP")]=eij[("ASP","GLN")]=-0.25
eij[("GLN","GLU")]=eij[("GLU","GLN")]=-0.12
eij[("GLN","HIS")]=eij[("HIS","GLN")]=0.06
eij[("GLN","ARG")]=eij[("ARG","GLN")]=0.04
eij[("GLN","LYS")]=eij[("LYS","GLN")]=0.18
eij[("GLN","PRO")]=eij[("PRO","GLN")]=0.11
eij[("ASP","ASP")]=0.21
eij[("ASP","GLU")]=eij[("GLU","ASP")]=0.68
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-0.53
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-0.26
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-0.09
eij[("ASP","PRO")]=eij[("PRO","ASP")]=0.33
eij[("GLU","GLU")]=0.6
eij[("GLU","HIS")]=eij[("HIS","GLU")]=-0.06
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-0.15
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-0.09
eij[("GLU","PRO")]=eij[("PRO","GLU")]=0.84
eij[("HIS","HIS")]=0.14
eij[("HIS","ARG")]=eij[("ARG","HIS")]=-0.01
eij[("HIS","LYS")]=eij[("LYS","HIS")]=0.14
eij[("HIS","PRO")]=eij[("PRO","HIS")]=-0.22
eij[("ARG","ARG")]=0.23
eij[("ARG","LYS")]=eij[("LYS","ARG")]=0.3
eij[("ARG","PRO")]=eij[("PRO","ARG")]=-0.02
eij[("LYS","LYS")]=1.45
eij[("LYS","PRO")]=eij[("PRO","LYS")]=0.51
eij[("PRO","PRO")]=0.28

def scoreEE(pdbid, ctype):

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
        
        #Eij=eij+err-eir-ejr	
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
        
    print("scoreEE:",str(round(Etot,2))," ctype:",ctype)
    return round(Etot,2)
