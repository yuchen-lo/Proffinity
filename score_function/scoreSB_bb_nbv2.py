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
eij[("CYS","CYS")]=-5.44
eij[("CYS","MET")]=eij[("MET","CYS")]=-4.99
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-5.8
eij[("CYS","ILE")]=eij[("ILE","CYS")]=-5.5
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-5.83
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-4.96
eij[("CYS","TRP")]=eij[("TRP","CYS")]=-4.95
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-4.16
eij[("CYS","ALA")]=eij[("ALA","CYS")]=-3.57
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-3.16
eij[("CYS","THR")]=eij[("THR","CYS")]=-3.11
eij[("CYS","SER")]=eij[("SER","CYS")]=-2.86
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-2.59
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-2.85
eij[("CYS","ASP")]=eij[("ASP","CYS")]=-2.41
eij[("CYS","GLU")]=eij[("GLU","CYS")]=-2.27
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-3.6
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-2.57
eij[("CYS","LYS")]=eij[("LYS","CYS")]=-1.95
eij[("CYS","PRO")]=eij[("PRO","CYS")]=-3.07
eij[("MET","MET")]=-5.46
eij[("MET","PHE")]=eij[("PHE","MET")]=-6.56
eij[("MET","ILE")]=eij[("ILE","MET")]=-6.02
eij[("MET","LEU")]=eij[("LEU","MET")]=-6.41
eij[("MET","VAL")]=eij[("VAL","MET")]=-5.32
eij[("MET","TRP")]=eij[("TRP","MET")]=-5.55
eij[("MET","TYR")]=eij[("TYR","MET")]=-4.91
eij[("MET","ALA")]=eij[("ALA","MET")]=-3.94
eij[("MET","GLY")]=eij[("GLY","MET")]=-3.39
eij[("MET","THR")]=eij[("THR","MET")]=-3.51
eij[("MET","SER")]=eij[("SER","MET")]=-3.03
eij[("MET","ASN")]=eij[("ASN","MET")]=-2.95
eij[("MET","GLN")]=eij[("GLN","MET")]=-3.30
eij[("MET","ASP")]=eij[("ASP","MET")]=-2.57
eij[("MET","GLU")]=eij[("GLU","MET")]=-2.89
eij[("MET","HIS")]=eij[("HIS","MET")]=-3.98
eij[("MET","ARG")]=eij[("ARG","MET")]=-3.12
eij[("MET","LYS")]=eij[("LYS","MET")]=-2.48
eij[("MET","PRO")]=eij[("PRO","MET")]=-3.45
eij[("PHE","PHE")]=-7.26
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-6.84
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-7.28
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-6.29
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-6.16
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-5.66
eij[("PHE","ALA")]=eij[("ALA","PHE")]=-4.81
eij[("PHE","GLY")]=eij[("GLY","PHE")]=-4.13
eij[("PHE","THR")]=eij[("THR","PHE")]=-4.28
eij[("PHE","SER")]=eij[("SER","PHE")]=-4.02
eij[("PHE","ASN")]=eij[("ASN","PHE")]=-3.75
eij[("PHE","GLN")]=eij[("GLN","PHE")]=-4.10
eij[("PHE","ASP")]=eij[("ASP","PHE")]=-3.48
eij[("PHE","GLU")]=eij[("GLU","PHE")]=-3.56
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-4.77
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-3.98
eij[("PHE","LYS")]=eij[("LYS","PHE")]=-3.36
eij[("PHE","PRO")]=eij[("PRO","PHE")]=-4.25
eij[("ILE","ILE")]=-6.54
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-7.04
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-6.05
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-5.78
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-5.25
eij[("ILE","ALA")]=eij[("ALA","ILE")]=-4.58
eij[("ILE","GLY")]=eij[("GLY","ILE")]=-3.78
eij[("ILE","THR")]=eij[("THR","ILE")]=-4.03
eij[("ILE","SER")]=eij[("SER","ILE")]=-3.52
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-3.24
eij[("ILE","GLN")]=eij[("GLN","ILE")]=-3.67
eij[("ILE","ASP")]=eij[("ASP","ILE")]=-3.17
eij[("ILE","GLU")]=eij[("GLU","ILE")]=-3.27
eij[("ILE","HIS")]=eij[("HIS","ILE")]=-4.14
eij[("ILE","ARG")]=eij[("ARG","ILE")]=-3.63
eij[("ILE","LYS")]=eij[("LYS","ILE")]=-3.01
eij[("ILE","PRO")]=eij[("PRO","ILE")]=-3.76
eij[("LEU","LEU")]=-7.37
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-6.48
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-6.14
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-5.67
eij[("LEU","ALA")]=eij[("ALA","LEU")]=-4.91
eij[("LEU","GLY")]=eij[("GLY","LEU")]=-4.16
eij[("LEU","THR")]=eij[("THR","LEU")]=-4.34
eij[("LEU","SER")]=eij[("SER","LEU")]=-3.92
eij[("LEU","ASN")]=eij[("ASN","LEU")]=-3.74
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-4.04
eij[("LEU","ASP")]=eij[("ASP","LEU")]=-3.40
eij[("LEU","GLU")]=eij[("GLU","LEU")]=-3.59
eij[("LEU","HIS")]=eij[("HIS","LEU")]=-4.54
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-4.03
eij[("LEU","LYS")]=eij[("LYS","LEU")]=-3.37
eij[("LEU","PRO")]=eij[("PRO","LEU")]=-4.20
eij[("VAL","VAL")]=-5.52
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-5.18
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-4.62
eij[("VAL","ALA")]=eij[("ALA","VAL")]=-4.04
eij[("VAL","GLY")]=eij[("GLY","VAL")]=-3.38
eij[("VAL","THR")]=eij[("THR","VAL")]=-3.46
eij[("VAL","SER")]=eij[("SER","VAL")]=-3.05
eij[("VAL","ASN")]=eij[("ASN","VAL")]=-2.83
eij[("VAL","GLN")]=eij[("GLN","VAL")]=-3.07
eij[("VAL","ASP")]=eij[("ASP","VAL")]=-2.48
eij[("VAL","GLU")]=eij[("GLU","VAL")]=-2.67
eij[("VAL","HIS")]=eij[("HIS","VAL")]=-3.58
eij[("VAL","ARG")]=eij[("ARG","VAL")]=-3.07
eij[("VAL","LYS")]=eij[("LYS","VAL")]=-2.49
eij[("VAL","PRO")]=eij[("PRO","VAL")]=-3.32
eij[("TRP","TRP")]=-5.06
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-4.66
eij[("TRP","ALA")]=eij[("ALA","TRP")]=-3.82
eij[("TRP","GLY")]=eij[("GLY","TRP")]=-3.42
eij[("TRP","THR")]=eij[("THR","TRP")]=-3.22
eij[("TRP","SER")]=eij[("SER","TRP")]=-2.99
eij[("TRP","ASN")]=eij[("ASN","TRP")]=-3.07
eij[("TRP","GLN")]=eij[("GLN","TRP")]=-3.11
eij[("TRP","ASP")]=eij[("ASP","TRP")]=-2.84
eij[("TRP","GLU")]=eij[("GLU","TRP")]=-2.99
eij[("TRP","HIS")]=eij[("HIS","TRP")]=-3.98
eij[("TRP","ARG")]=eij[("ARG","TRP")]=-3.41
eij[("TRP","LYS")]=eij[("LYS","TRP")]=-2.69
eij[("TRP","PRO")]=eij[("PRO","TRP")]=-3.73
eij[("TYR","TYR")]=-4.17
eij[("TYR","ALA")]=eij[("ALA","TYR")]=-3.36
eij[("TYR","GLY")]=eij[("GLY","TYR")]=-3.01
eij[("TYR","THR")]=eij[("THR","TYR")]=-3.01
eij[("TYR","SER")]=eij[("SER","TYR")]=-2.78
eij[("TYR","ASN")]=eij[("ASN","TYR")]=-2.76
eij[("TYR","GLN")]=eij[("GLN","TYR")]=-2.97
eij[("TYR","ASP")]=eij[("ASP","TYR")]=-2.76
eij[("TYR","GLU")]=eij[("GLU","TYR")]=-2.79
eij[("TYR","HIS")]=eij[("HIS","TYR")]=-3.52
eij[("TYR","ARG")]=eij[("ARG","TYR")]=-3.16
eij[("TYR","LYS")]=eij[("LYS","TYR")]=-2.60
eij[("TYR","PRO")]=eij[("PRO","TYR")]=-3.19
eij[("ALA","ALA")]=-2.72
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-2.31
eij[("ALA","THR")]=eij[("THR","ALA")]=-2.32
eij[("ALA","SER")]=eij[("SER","ALA")]=-2.01
eij[("ALA","ASN")]=eij[("ASN","ALA")]=-1.84
eij[("ALA","GLN")]=eij[("GLN","ALA")]=-1.89
eij[("ALA","ASP")]=eij[("ASP","ALA")]=-1.70
eij[("ALA","GLU")]=eij[("GLU","ALA")]=-1.51
eij[("ALA","HIS")]=eij[("HIS","ALA")]=-2.41
eij[("ALA","ARG")]=eij[("ARG","ALA")]=-1.83
eij[("ALA","LYS")]=eij[("LYS","ALA")]=-1.31
eij[("ALA","PRO")]=eij[("PRO","ALA")]=-2.03
eij[("GLY","GLY")]=-2.24
eij[("GLY","THR")]=eij[("THR","GLY")]=-2.08
eij[("GLY","SER")]=eij[("SER","GLY")]=-1.82
eij[("GLY","ASN")]=eij[("ASN","GLY")]=-1.74
eij[("GLY","GLN")]=eij[("GLN","GLY")]=-1.66
eij[("GLY","ASP")]=eij[("ASP","GLY")]=-1.59
eij[("GLY","GLU")]=eij[("GLU","GLY")]=-1.22
eij[("GLY","HIS")]=eij[("HIS","GLY")]=-2.15
eij[("GLY","ARG")]=eij[("ARG","GLY")]=-1.72
eij[("GLY","LYS")]=eij[("LYS","GLY")]=-1.15
eij[("GLY","PRO")]=eij[("PRO","GLY")]=-1.87
eij[("THR","THR")]=-2.12
eij[("THR","SER")]=eij[("SER","THR")]=-1.96
eij[("THR","ASN")]=eij[("ASN","THR")]=-1.88
eij[("THR","GLN")]=eij[("GLN","THR")]=-1.90
eij[("THR","ASP")]=eij[("ASP","THR")]=-1.80
eij[("THR","GLU")]=eij[("GLU","THR")]=-1.74
eij[("THR","HIS")]=eij[("HIS","THR")]=-2.42
eij[("THR","ARG")]=eij[("ARG","THR")]=-1.90
eij[("THR","LYS")]=eij[("LYS","THR")]=-1.31
eij[("THR","PRO")]=eij[("PRO","THR")]=-1.90
eij[("SER","SER")]=-1.67
eij[("SER","ASN")]=eij[("ASN","SER")]=-1.58
eij[("SER","GLN")]=eij[("GLN","SER")]=-1.49
eij[("SER","ASP")]=eij[("ASP","SER")]=-1.63
eij[("SER","GLU")]=eij[("GLU","SER")]=-1.48
eij[("SER","HIS")]=eij[("HIS","SER")]=-2.11
eij[("SER","ARG")]=eij[("ARG","SER")]=-1.62
eij[("SER","LYS")]=eij[("LYS","SER")]=-1.05
eij[("SER","PRO")]=eij[("PRO","SER")]=-1.57
eij[("ASN","ASN")]=-1.68
eij[("ASN","GLN")]=eij[("GLN","ASN")]=-1.71
eij[("ASN","ASP")]=eij[("ASP","ASN")]=-1.68
eij[("ASN","GLU")]=eij[("GLU","ASN")]=-1.51
eij[("ASN","HIS")]=eij[("HIS","ASN")]=-2.08
eij[("ASN","ARG")]=eij[("ARG","ASN")]=-1.64
eij[("ASN","LYS")]=eij[("LYS","ASN")]=-1.21
eij[("ASN","PRO")]=eij[("PRO","ASN")]=-1.53
eij[("GLN","GLN")]=-1.54
eij[("GLN","ASP")]=eij[("ASP","GLN")]=-1.46
eij[("GLN","GLU")]=eij[("GLU","GLN")]=-1.42
eij[("GLN","HIS")]=eij[("HIS","GLN")]=-1.98
eij[("GLN","ARG")]=eij[("ARG","GLN")]=-1.80
eij[("GLN","LYS")]=eij[("LYS","GLN")]=-1.29
eij[("GLN","PRO")]=eij[("PRO","GLN")]=-1.73
eij[("ASP","ASP")]=-1.21
eij[("ASP","GLU")]=eij[("GLU","ASP")]=-1.02
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-2.32
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-2.29
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-1.68
eij[("ASP","PRO")]=eij[("PRO","ASP")]=-1.33
eij[("GLU","GLU")]=-0.91
eij[("GLU","HIS")]=eij[("HIS","GLU")]=-2.15
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-2.27
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-1.8
eij[("GLU","PRO")]=eij[("PRO","GLU")]=-1.26
eij[("HIS","HIS")]=-3.05
eij[("HIS","ARG")]=eij[("ARG","HIS")]=-2.16
eij[("HIS","LYS")]=eij[("LYS","HIS")]=-1.35
eij[("HIS","PRO")]=eij[("PRO","HIS")]=-2.25
eij[("ARG","ARG")]=-1.55
eij[("ARG","LYS")]=eij[("LYS","ARG")]=-0.59
eij[("ARG","PRO")]=eij[("PRO","ARG")]=-1.70
eij[("LYS","LYS")]=-0.12
eij[("LYS","PRO")]=eij[("PRO","LYS")]=-0.97
eij[("PRO","PRO")]=-1.75

err=-2.55


ekr={}
ekr['CYS']=-3.57
ekr['MET']=-3.92
ekr['PHE']=-4.76
ekr['ILE']=-4.42
ekr['LEU']=-4.81
ekr['VAL']=-3.89
ekr['TRP']=-3.81
ekr['TYR']=-3.41
ekr['ALA']=-2.57
ekr['GLY']=-2.19
ekr['THR']=-2.29
ekr['SER']=-1.98
ekr['ASN']=-1.92
ekr['GLN']=-2
ekr['ASP']=-1.84
ekr['GLU']=-1.79
ekr['HIS']=-2.56
ekr['ARG']=-2.11
ekr['LYS']=-1.52
ekr['PRO']=-2.09


#Repusive Energy
Q={}
Q['CYS']=6.65
Q['MET']=6.14
Q['PHE']=5.87
Q['ILE']=6.04
Q['LEU']=6.09
Q['VAL']=6.16
Q['TRP']=5.79
Q['TYR']=6.04
Q['ALA']=6.33
Q['GLY']=6.28
Q['THR']=6.49
Q['SER']=6.58
Q['ASN']=6.57
Q['GLN']=6.47
Q['ASP']=6.49
Q['GLU']=6.24
Q['HIS']=6.24
Q['ARG']=6.32
Q['LYS']=6.57
Q['PRO']=5.86

#Distribution for the total number of residues type surrounded by N residues
def N(resi,n):
    ntot=0
    if resi == 'CYS':
        ntot= -0.85*n**4 + 31.94*n**3 - 432.94*n**2 + 2445.84*n - 4655.44
    elif resi == 'MET':
        ntot= -0.70*n**4 + 27.41*n**3 - 385.76*n**2 + 2272.82*n -4542.80
    elif resi == 'PHE':
        ntot= -1.34*n**4 + 51.37*n**3 - 708.65*n**2 + 4064.01*n - 7774.64
    elif resi == 'ILE':
        ntot= -1.93*n**4 + 74.99*n**3 - 1049.10*n**2 + 6127.05*n - 12054.42
    elif resi == 'LEU':
        ntot= -3.14*n**4 + 123.89*n**3 - 1762.22*n**2 + 10508.19*n - 21302.37
    elif resi == 'VAL':
        ntot= -2.71*n**4 + 104.64*n**3 - 1456.43*n**2 + 8479.84*n - 16725.74
    elif resi == 'TRP':
        ntot= -0.38*n**4 + 13.94*n**3 - 180.72*n**2 + 945.15*n - 1523.79
    elif resi == 'TYR':
        ntot= -1.02*n**4 + 37.71*n**3 - 495.50*n**2 + 2657.41*n - 4545.06
    elif resi == 'ALA':
        ntot= -2.4*n**4 + 89.11*n**3 - 1179.69*n**2 + 6447.48*n - 11611.99
    elif resi == 'GLY':
        ntot= -0.69*n**4 + 23.63*n**3 - 273.59*n**2 + 1121.55*n - 608.88
    elif resi == 'THR':
        ntot= -0.46*n**4 + 16.31*n**3 - 197.86*n**2 + 877.41*n - 706.10
    elif resi == 'SER':
        ntot= -0.77*n**4 + 27.53*n**3 - 345.14*n**2 + 1701.3*n - 2365.65
    elif resi == 'ASN':
        ntot= -0.09*n**4 + 2.15*n**3 - 3.39*n**2 - 223.97*n + 1324.07
    elif resi == 'GLN':
        ntot= -0.14*n**4 + 4.49*n**3 - 45.2*n**2 + 102.37*n + 388.83
    elif resi == 'ASP':
        ntot= 0.054*n**4 - 2.97*n**3 + 62.6*n**2 - 591.25*n + 2087.63
    elif resi == 'GLU':
        ntot= -0.177*n**4 + 5.45*n**3 - 49.44*n**2 + 57.34*n + 687.12
    elif resi == 'HIS':
        ntot= -0.18*n**4 + 6.39*n**3 - 76.94*n**2 + 332.7*n - 222.09
    elif resi == 'ARG':
        ntot= -0.13*n**4 + 3.16*n**3 - 11.96*n**2 - 207.58*n + 1362.7
    elif resi == 'LYS':
        ntot= 0.16*n**4 - 7.48*n**3 + 126.52*n**2 - 949.43*n + 2667.44
    elif resi == 'PRO':
        ntot= -0.14*n**4 + 3.99*n**3 - 30.28*n**2 - 52.98*n + 930.03
        
    if ntot <0:
        ntot=0    
        
    return ntot
      

def scoreSB_BB_NB(pdbid, ctype):
# Define the parser

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
    rpack = 0
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

        Eij=-1*(eij[(resA,resB)]+err-ekr[resA]-ekr[resB])
    
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
            rpack=rpack+1
        else:
            last_line=0   
	
        if (resAid != resPid or last_line==1):
            if (rpack >= Q[resA]):
                Ep=(Q[resA]/rpack-1)*Ec-log((N(resA,rpack)+10**-6)/(N(resA,Q[resA])+10**-6))
            else:
                Ep=0
            
            Etot=Etot+Ec+Ep
            #print("resPid:"+resPid+" Etot:"+str(round(Etot,2))+" Ec:"+str(round(Ec,2))+" lines:"+str(num_lines)+"/"+str(total_lines))
            Ep=0
            Ec=0
            rpack=0
            resPid=resAid

        Ec=Ec+Eij
        rpack=rpack+1


    print("scoreSB_BB_NB:",str(round(Etot/100,2))," ctype:",ctype)
    return round(Etot/100,2)
