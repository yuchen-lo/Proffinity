import argparse
import sys
from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt, log, exp

#################
#Xu et al. (1998)
#e1=e1+err-eir-ejr
e1={}
e1[("ALA","ALA")]=-14
e1[("ALA","ARG")]=e1[("ARG","ALA")]=27
e1[("ALA","ASN")]=e1[("ASN","ALA")]=10
e1[("ALA","ASP")]=e1[("ASP","ALA")]=22
e1[("ALA","CYS")]=e1[("CYS","ALA")]=33
e1[("ALA","GLN")]=e1[("GLN","ALA")]=3
e1[("ALA","GLU")]=e1[("GLU","ALA")]=12
e1[("ALA","GLY")]=e1[("GLY","ALA")]=1
e1[("ALA","HIS")]=e1[("HIS","ALA")]=6
e1[("ALA","ILE")]=e1[("ILE","ALA")]=-11
e1[("ALA","LEU")]=e1[("LEU","ALA")]=-18
e1[("ALA","LYS")]=e1[("LYS","ALA")]=12
e1[("ALA","MET")]=e1[("MET","ALA")]=-7
e1[("ALA","PHE")]=e1[("PHE","ALA")]=-6
e1[("ALA","PRO")]=e1[("PRO","ALA")]=17
e1[("ALA","SER")]=e1[("SER","ALA")]=17
e1[("ALA","THR")]=e1[("THR","ALA")]=6
e1[("ALA","TRP")]=e1[("TRP","ALA")]=5
e1[("ALA","TYR")]=e1[("TYR","ALA")]=5
e1[("ALA","VAL")]=e1[("VAL","ALA")]=-11
e1[("ARG","ARG")]=-2
e1[("ARG","ASN")]=e1[("ASN","ARG")]=-9
e1[("ARG","ASP")]=e1[("ASP","ARG")]=-62
e1[("ARG","CYS")]=e1[("CYS","ARG")]=7
e1[("ARG","GLN")]=e1[("GLN","ARG")]=-6
e1[("ARG","GLU")]=e1[("GLU","ARG")]=-56
e1[("ARG","GLY")]=e1[("GLY","ARG")]=-8
e1[("ARG","HIS")]=e1[("HIS","ARG")]=-26
e1[("ARG","ILE")]=e1[("ILE","ARG")]=11
e1[("ARG","LEU")]=e1[("LEU","ARG")]=26
e1[("ARG","LYS")]=e1[("LYS","ARG")]=31
e1[("ARG","MET")]=e1[("MET","ARG")]=30
e1[("ARG","PHE")]=e1[("PHE","ARG")]=6
e1[("ARG","PRO")]=e1[("PRO","ARG")]=-3
e1[("ARG","SER")]=e1[("SER","ARG")]=-8
e1[("ARG","THR")]=e1[("THR","ARG")]=6
e1[("ARG","TRP")]=e1[("TRP","ARG")]=-15
e1[("ARG","TYR")]=e1[("TYR","ARG")]=-13
e1[("ARG","VAL")]=e1[("VAL","ARG")]=17
e1[("ASN","ASN")]=-43
e1[("ASN","ASP")]=e1[("ASP","ASN")]=-42
e1[("ASN","CYS")]=e1[("CYS","ASN")]=11
e1[("ASN","GLN")]=e1[("GLN","ASN")]=-20
e1[("ASN","GLU")]=e1[("GLU","ASN")]=-14
e1[("ASN","GLY")]=e1[("GLY","ASN")]=-10
e1[("ASN","HIS")]=e1[("HIS","ASN")]=6
e1[("ASN","ILE")]=e1[("ILE","ASN")]=35
e1[("ASN","LEU")]=e1[("LEU","ASN")]=36
e1[("ASN","LYS")]=e1[("LYS","ASN")]=-20
e1[("ASN","MET")]=e1[("MET","ASN")]=31
e1[("ASN","PHE")]=e1[("PHE","ASN")]=20
e1[("ASN","PRO")]=e1[("PRO","ASN")]=-21
e1[("ASN","SER")]=e1[("SER","ASN")]=-22
e1[("ASN","THR")]=e1[("THR","ASN")]=-23
e1[("ASN","TRP")]=e1[("TRP","ASN")]=-2
e1[("ASN","TYR")]=e1[("TYR","ASN")]=5
e1[("ASN","VAL")]=e1[("VAL","ASN")]=30
e1[("ASP","ASP")]=2
e1[("ASP","CYS")]=e1[("CYS","ASP")]=28
e1[("ASP","GLN")]=e1[("GLN","ASP")]=7
e1[("ASP","GLU")]=e1[("GLU","ASP")]=14
e1[("ASP","GLY")]=e1[("GLY","ASP")]=-27
e1[("ASP","HIS")]=e1[("HIS","ASP")]=-45
e1[("ASP","ILE")]=e1[("ILE","ASP")]=32
e1[("ASP","LEU")]=e1[("LEU","ASP")]=37
e1[("ASP","LYS")]=e1[("LYS","ASP")]=-56
e1[("ASP","MET")]=e1[("MET","ASP")]=21
e1[("ASP","PHE")]=e1[("PHE","ASP")]=28
e1[("ASP","PRO")]=e1[("PRO","ASP")]=-3
e1[("ASP","SER")]=e1[("SER","ASP")]=-30
e1[("ASP","THR")]=e1[("THR","ASP")]=-20
e1[("ASP","TRP")]=e1[("TRP","ASP")]=10
e1[("ASP","TYR")]=e1[("TYR","ASP")]=27
e1[("ASP","VAL")]=e1[("VAL","ASP")]=43
e1[("CYS","CYS")]=-192
e1[("CYS","GLN")]=e1[("GLN","CYS")]=19
e1[("CYS","GLU")]=e1[("GLU","CYS")]=12
e1[("CYS","GLY")]=e1[("GLY","CYS")]=9
e1[("CYS","HIS")]=e1[("HIS","CYS")]=19
e1[("CYS","ILE")]=e1[("ILE","CYS")]=15
e1[("CYS","LEU")]=e1[("LEU","CYS")]=24
e1[("CYS","LYS")]=e1[("LYS","CYS")]=25
e1[("CYS","MET")]=e1[("MET","CYS")]=5
e1[("CYS","PHE")]=e1[("PHE","CYS")]=3
e1[("CYS","PRO")]=e1[("PRO","CYS")]=11
e1[("CYS","SER")]=e1[("SER","CYS")]=1
e1[("CYS","THR")]=e1[("THR","CYS")]=37
e1[("CYS","TRP")]=e1[("TRP","CYS")]=5
e1[("CYS","TYR")]=e1[("TYR","CYS")]=6
e1[("CYS","VAL")]=e1[("VAL","CYS")]=20
e1[("GLN","GLN")]=-11
e1[("GLN","GLU")]=e1[("GLU","GLN")]=1
e1[("GLN","GLY")]=e1[("GLY","GLN")]=-7
e1[("GLN","HIS")]=e1[("HIS","GLN")]=27
e1[("GLN","ILE")]=e1[("ILE","GLN")]=24
e1[("GLN","LEU")]=e1[("LEU","GLN")]=2
e1[("GLN","LYS")]=e1[("LYS","GLN")]=-18
e1[("GLN","MET")]=e1[("MET","GLN")]=3
e1[("GLN","PHE")]=e1[("PHE","GLN")]=7
e1[("GLN","PRO")]=e1[("PRO","GLN")]=-8
e1[("GLN","SER")]=e1[("SER","GLN")]=-16
e1[("GLN","THR")]=e1[("THR","GLN")]=-15
e1[("GLN","TRP")]=e1[("TRP","GLN")]=-1
e1[("GLN","TYR")]=e1[("TYR","GLN")]=-9
e1[("GLN","VAL")]=e1[("VAL","GLN")]=18
e1[("GLU","GLU")]=7
e1[("GLU","GLY")]=e1[("GLY","GLU")]=-3
e1[("GLU","HIS")]=e1[("HIS","GLU")]=-37
e1[("GLU","ILE")]=e1[("ILE","GLU")]=29
e1[("GLU","LEU")]=e1[("LEU","GLU")]=25
e1[("GLU","LYS")]=e1[("LYS","GLU")]=-67
e1[("GLU","MET")]=e1[("MET","GLU")]=14
e1[("GLU","PHE")]=e1[("PHE","GLU")]=23
e1[("GLU","PRO")]=e1[("PRO","GLU")]=-10
e1[("GLU","SER")]=e1[("SER","GLU")]=-21
e1[("GLU","THR")]=e1[("THR","GLU")]=-21
e1[("GLU","TRP")]=e1[("TRP","GLU")]=16
e1[("GLU","TYR")]=e1[("TYR","GLU")]=27
e1[("GLU","VAL")]=e1[("VAL","GLU")]=23
e1[("GLY","GLY")]=-29
e1[("GLY","HIS")]=e1[("HIS","GLY")]=7
e1[("GLY","ILE")]=e1[("ILE","GLY")]=18
e1[("GLY","LEU")]=e1[("LEU","GLY")]=24
e1[("GLY","LYS")]=e1[("LYS","GLY")]=10
e1[("GLY","MET")]=e1[("MET","GLY")]=1
e1[("GLY","PHE")]=e1[("PHE","GLY")]=11
e1[("GLY","PRO")]=e1[("PRO","GLY")]=-7
e1[("GLY","SER")]=e1[("SER","GLY")]=-19
e1[("GLY","THR")]=e1[("THR","GLY")]=-7
e1[("GLY","TRP")]=e1[("TRP","GLY")]=-7
e1[("GLY","TYR")]=e1[("TYR","GLY")]=6
e1[("GLY","VAL")]=e1[("VAL","GLY")]=20
e1[("HIS","HIS")]=-45
e1[("HIS","ILE")]=e1[("ILE","HIS")]=29
e1[("HIS","LEU")]=e1[("LEU","HIS")]=20
e1[("HIS","LYS")]=e1[("LYS","HIS")]=5
e1[("HIS","MET")]=e1[("MET","HIS")]=-1
e1[("HIS","PHE")]=e1[("PHE","HIS")]=16
e1[("HIS","PRO")]=e1[("PRO","HIS")]=-6
e1[("HIS","SER")]=e1[("SER","HIS")]=-13
e1[("HIS","THR")]=e1[("THR","HIS")]=-24
e1[("HIS","TRP")]=e1[("TRP","HIS")]=-21
e1[("HIS","TYR")]=e1[("TYR","HIS")]=3
e1[("HIS","VAL")]=e1[("VAL","HIS")]=20
e1[("ILE","ILE")]=-33
e1[("ILE","LEU")]=e1[("LEU","ILE")]=-16
e1[("ILE","LYS")]=e1[("LYS","ILE")]=19
e1[("ILE","MET")]=e1[("MET","ILE")]=-1
e1[("ILE","PHE")]=e1[("PHE","ILE")]=-10
e1[("ILE","PRO")]=e1[("PRO","ILE")]=37
e1[("ILE","SER")]=e1[("SER","ILE")]=21
e1[("ILE","THR")]=e1[("THR","ILE")]=11
e1[("ILE","TRP")]=e1[("TRP","ILE")]=-2
e1[("ILE","TYR")]=e1[("TYR","ILE")]=-16
e1[("ILE","VAL")]=e1[("VAL","ILE")]=-23
e1[("LEU","LEU")]=-28
e1[("LEU","LYS")]=e1[("LYS","LEU")]=18
e1[("LEU","MET")]=e1[("MET","LEU")]=-11
e1[("LEU","PHE")]=e1[("PHE","LEU")]=-20
e1[("LEU","PRO")]=e1[("PRO","LEU")]=22
e1[("LEU","SER")]=e1[("SER","LEU")]=27
e1[("LEU","THR")]=e1[("THR","LEU")]=22
e1[("LEU","TRP")]=e1[("TRP","LEU")]=8
e1[("LEU","TYR")]=e1[("TYR","LEU")]=-9
e1[("LEU","VAL")]=e1[("VAL","LEU")]=-22
e1[("LYS","LYS")]=12
e1[("LYS","MET")]=e1[("MET","LYS")]=30
e1[("LYS","PHE")]=e1[("PHE","LYS")]=-2
e1[("LYS","PRO")]=e1[("PRO","LYS")]=-5
e1[("LYS","SER")]=e1[("SER","LYS")]=-6
e1[("LYS","THR")]=e1[("THR","LYS")]=-2
e1[("LYS","TRP")]=e1[("TRP","LYS")]=3
e1[("LYS","TYR")]=e1[("TYR","LYS")]=-31
e1[("LYS","VAL")]=e1[("VAL","LYS")]=27
e1[("MET","MET")]=-49
e1[("MET","PHE")]=e1[("PHE","MET")]=-27
e1[("MET","PRO")]=e1[("PRO","MET")]=3
e1[("MET","SER")]=e1[("SER","MET")]=19
e1[("MET","THR")]=e1[("THR","MET")]=16
e1[("MET","TRP")]=e1[("TRP","MET")]=-1
e1[("MET","TYR")]=e1[("TYR","MET")]=-17
e1[("MET","VAL")]=e1[("VAL","MET")]=-5
e1[("PHE","PHE")]=-21
e1[("PHE","PRO")]=e1[("PRO","PHE")]=-2
e1[("PHE","SER")]=e1[("SER","PHE")]=11
e1[("PHE","THR")]=e1[("THR","PHE")]=28
e1[("PHE","TRP")]=e1[("TRP","PHE")]=3
e1[("PHE","TYR")]=e1[("TYR","PHE")]=-1
e1[("PHE","VAL")]=e1[("VAL","PHE")]=-4
e1[("PRO","PRO")]=-21
e1[("PRO","SER")]=e1[("SER","PRO")]=-16
e1[("PRO","THR")]=e1[("THR","PRO")]=-10
e1[("PRO","TRP")]=e1[("TRP","PRO")]=-43
e1[("PRO","TYR")]=e1[("TYR","PRO")]=-8
e1[("PRO","VAL")]=e1[("VAL","PRO")]=5
e1[("SER","SER")]=-18
e1[("SER","THR")]=e1[("THR","SER")]=-22
e1[("SER","TRP")]=e1[("TRP","SER")]=13
e1[("SER","TYR")]=e1[("TYR","SER")]=10
e1[("SER","VAL")]=e1[("VAL","SER")]=27
e1[("THR","THR")]=-21
e1[("THR","TRP")]=e1[("TRP","THR")]=9
e1[("THR","TYR")]=e1[("TYR","THR")]=16
e1[("THR","VAL")]=e1[("VAL","THR")]=7
e1[("TRP","TRP")]=-2
e1[("TRP","TYR")]=e1[("TYR","TRP")]=-10
e1[("TRP","VAL")]=e1[("VAL","TRP")]=10
e1[("TYR","TYR")]=-1
e1[("TYR","VAL")]=e1[("VAL","TYR")]=11
e1[("VAL","VAL")]=-32


def scoreKP(pdbid, ctype):

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


    resPid = Lines[0].split(',')[2]
    Ec=0
    Etot=0
    num_lines=0
    last_line=0

    # Strips the newline character
    for line in Lines:
        resA=line.split(',')[1]
        resAid=line.split(',')[2]
        resB=line.split(',')[4]
        dist=float(line.split(',')[6])
        num_lines=num_lines+1
          
        e1k=e1[(resA,resB)]
        Eij=-1*(e1k/100)
    
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

    print("scoreKP:",str(round(Etot,2))," ctype:",ctype)
    return round(Etot,2)
