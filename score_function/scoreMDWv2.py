import argparse
import sys
from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt, log

#################
#Kamagata et al. (2022)
#Eij=eij+err-eir-ejr
eij={}
eij[("CYS","CYS")]=-1.92478
eij[("CYS","MET")]=eij[("MET","CYS")]=-1.66434
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-1.69439
eij[("CYS","ILE")]=eij[("ILE","CYS")]=-1.79699
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-1.32278
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-2.07097
eij[("CYS","TRP")]=eij[("TRP","CYS")]=-2.87560
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-1.90875
eij[("CYS","ALA")]=eij[("ALA","CYS")]=-2.12923
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-2.00976
eij[("CYS","THR")]=eij[("THR","CYS")]=-2.28725
eij[("CYS","SER")]=eij[("SER","CYS")]=-2.17174
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-2.32682
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-1.48432
eij[("CYS","ASP")]=eij[("ASP","CYS")]=-1.81920
eij[("CYS","GLU")]=eij[("GLU","CYS")]=-1.58849
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-1.30303
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-1.78557
eij[("CYS","LYS")]=eij[("LYS","CYS")]=-1.71650
eij[("CYS","PRO")]=eij[("PRO","CYS")]=-1.50752
eij[("MET","MET")]=-1.58835
eij[("MET","PHE")]=eij[("PHE","MET")]=-1.41570
eij[("MET","ILE")]=eij[("ILE","MET")]=-1.55946
eij[("MET","LEU")]=eij[("LEU","MET")]=-1.33137
eij[("MET","VAL")]=eij[("VAL","MET")]=-1.98223
eij[("MET","TRP")]=eij[("TRP","MET")]=-3.10936
eij[("MET","TYR")]=eij[("TYR","MET")]=-1.93880
eij[("MET","ALA")]=eij[("ALA","MET")]=-1.56745
eij[("MET","GLY")]=eij[("GLY","MET")]=-1.38147
eij[("MET","THR")]=eij[("THR","MET")]=-1.49344
eij[("MET","SER")]=eij[("SER","MET")]=-2.12528
eij[("MET","ASN")]=eij[("ASN","MET")]=-2.27205
eij[("MET","GLN")]=eij[("GLN","MET")]=-1.38862
eij[("MET","ASP")]=eij[("ASP","MET")]=-1.38297
eij[("MET","GLU")]=eij[("GLU","MET")]=-1.57434
eij[("MET","HIS")]=eij[("HIS","MET")]=-1.55421
eij[("MET","ARG")]=eij[("ARG","MET")]=-2.25171
eij[("MET","LYS")]=eij[("LYS","MET")]=-1.26644
eij[("MET","PRO")]=eij[("PRO","MET")]=-1.40234
eij[("PHE","PHE")]=-1.98447
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-1.49885
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-1.89611
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-1.30168
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-3.38140
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-2.77317
eij[("PHE","ALA")]=eij[("ALA","PHE")]=-1.83994
eij[("PHE","GLY")]=eij[("GLY","PHE")]=-1.56865
eij[("PHE","THR")]=eij[("THR","PHE")]=-1.79604
eij[("PHE","SER")]=eij[("SER","PHE")]=-1.66355
eij[("PHE","ASN")]=eij[("ASN","PHE")]=-2.40241
eij[("PHE","GLN")]=eij[("GLN","PHE")]=-2.03858
eij[("PHE","ASP")]=eij[("ASP","PHE")]=-1.69530
eij[("PHE","GLU")]=eij[("GLU","PHE")]=-1.09783
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-2.05199
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-3.03696
eij[("PHE","LYS")]=eij[("LYS","PHE")]=-1.30816
eij[("PHE","PRO")]=eij[("PRO","PHE")]=-1.97895
eij[("ILE","ILE")]=-1.46660
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-1.43020
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-1.43962
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-2.10208
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-2.25773
eij[("ILE","ALA")]=eij[("ALA","ILE")]=-1.48696
eij[("ILE","GLY")]=eij[("GLY","ILE")]=-1.38644
eij[("ILE","THR")]=eij[("THR","ILE")]=-1.56220
eij[("ILE","SER")]=eij[("SER","ILE")]=-1.66610
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-2.55811
eij[("ILE","GLN")]=eij[("GLN","ILE")]=-2.01845
eij[("ILE","ASP")]=eij[("ASP","ILE")]=-1.73345
eij[("ILE","GLU")]=eij[("GLU","ILE")]=-1.59699
eij[("ILE","HIS")]=eij[("HIS","ILE")]=-1.62433
eij[("ILE","ARG")]=eij[("ARG","ILE")]=-1.94227
eij[("ILE","LYS")]=eij[("LYS","ILE")]=-1.58615
eij[("ILE","PRO")]=eij[("PRO","ILE")]=-1.52216
eij[("LEU","LEU")]=-0.93259
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-1.63149
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-1.88308
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-1.93806
eij[("LEU","ALA")]=eij[("ALA","LEU")]=-1.54687
eij[("LEU","GLY")]=eij[("GLY","LEU")]=-1.84977
eij[("LEU","THR")]=eij[("THR","LEU")]=-2.23461
eij[("LEU","SER")]=eij[("SER","LEU")]=-2.40301
eij[("LEU","ASN")]=eij[("ASN","LEU")]=-2.27150
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-2.07635
eij[("LEU","ASP")]=eij[("ASP","LEU")]=-1.98325
eij[("LEU","GLU")]=eij[("GLU","LEU")]=-1.06194
eij[("LEU","HIS")]=eij[("HIS","LEU")]=-1.63835
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-2.18864
eij[("LEU","LYS")]=eij[("LYS","LEU")]=-1.67308
eij[("LEU","PRO")]=eij[("PRO","LEU")]=-1.63923
eij[("VAL","VAL")]=-2.26125
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-1.80166
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-1.73793
eij[("VAL","ALA")]=eij[("ALA","VAL")]=-1.53896
eij[("VAL","GLY")]=eij[("GLY","VAL")]=-1.38517
eij[("VAL","THR")]=eij[("THR","VAL")]=-1.78412
eij[("VAL","SER")]=eij[("SER","VAL")]=-2.46446
eij[("VAL","ASN")]=eij[("ASN","VAL")]=-2.57377
eij[("VAL","GLN")]=eij[("GLN","VAL")]=-1.95542
eij[("VAL","ASP")]=eij[("ASP","VAL")]=-1.70581
eij[("VAL","GLU")]=eij[("GLU","VAL")]=-1.56714
eij[("VAL","HIS")]=eij[("HIS","VAL")]=-2.32664
eij[("VAL","ARG")]=eij[("ARG","VAL")]=-1.82820
eij[("VAL","LYS")]=eij[("LYS","VAL")]=-1.54127
eij[("VAL","PRO")]=eij[("PRO","VAL")]=-0.98615
eij[("TRP","TRP")]=-3.99217
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-3.25044
eij[("TRP","ALA")]=eij[("ALA","TRP")]=-1.22216
eij[("TRP","GLY")]=eij[("GLY","TRP")]=-1.64207
eij[("TRP","THR")]=eij[("THR","TRP")]=-2.18168
eij[("TRP","SER")]=eij[("SER","TRP")]=-2.16924
eij[("TRP","ASN")]=eij[("ASN","TRP")]=-2.34652
eij[("TRP","GLN")]=eij[("GLN","TRP")]=-3.04111
eij[("TRP","ASP")]=eij[("ASP","TRP")]=-1.37952
eij[("TRP","GLU")]=eij[("GLU","TRP")]=-1.30757
eij[("TRP","HIS")]=eij[("HIS","TRP")]=-2.83422
eij[("TRP","ARG")]=eij[("ARG","TRP")]=-6.16511
eij[("TRP","LYS")]=eij[("LYS","TRP")]=-3.43548
eij[("TRP","PRO")]=eij[("PRO","TRP")]=-2.15874
eij[("TYR","TYR")]=-3.20608
eij[("TYR","ALA")]=eij[("ALA","TYR")]=-1.60138
eij[("TYR","GLY")]=eij[("GLY","TYR")]=-1.26081
eij[("TYR","THR")]=eij[("THR","TYR")]=-1.85248
eij[("TYR","SER")]=eij[("SER","TYR")]=-1.77031
eij[("TYR","ASN")]=eij[("ASN","TYR")]=-1.88370
eij[("TYR","GLN")]=eij[("GLN","TYR")]=-1.91097
eij[("TYR","ASP")]=eij[("ASP","TYR")]=-1.67571
eij[("TYR","GLU")]=eij[("GLU","TYR")]=-1.53305
eij[("TYR","HIS")]=eij[("HIS","TYR")]=-2.15948
eij[("TYR","ARG")]=eij[("ARG","TYR")]=-4.74134
eij[("TYR","LYS")]=eij[("LYS","TYR")]=-1.54497
eij[("TYR","PRO")]=eij[("PRO","TYR")]=-2.07329
eij[("ALA","ALA")]=-1.29831
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-1.77891
eij[("ALA","THR")]=eij[("THR","ALA")]=-2.00108
eij[("ALA","SER")]=eij[("SER","ALA")]=-2.20592
eij[("ALA","ASN")]=eij[("ASN","ALA")]=-2.63429
eij[("ALA","GLN")]=eij[("GLN","ALA")]=-1.70793
eij[("ALA","ASP")]=eij[("ASP","ALA")]=-2.44411
eij[("ALA","GLU")]=eij[("GLU","ALA")]=-1.35144
eij[("ALA","HIS")]=eij[("HIS","ALA")]=-1.93409
eij[("ALA","ARG")]=eij[("ARG","ALA")]=-1.23631
eij[("ALA","LYS")]=eij[("LYS","ALA")]=-1.48560
eij[("ALA","PRO")]=eij[("PRO","ALA")]=-1.78633
eij[("GLY","GLY")]=-1.87863
eij[("GLY","THR")]=eij[("THR","GLY")]=-1.93278
eij[("GLY","SER")]=eij[("SER","GLY")]=-2.14923
eij[("GLY","ASN")]=eij[("ASN","GLY")]=-2.92082
eij[("GLY","GLN")]=eij[("GLN","GLY")]=-1.71393
eij[("GLY","ASP")]=eij[("ASP","GLY")]=-2.52797
eij[("GLY","GLU")]=eij[("GLU","GLY")]=-0.89530
eij[("GLY","HIS")]=eij[("HIS","GLY")]=-1.80133
eij[("GLY","ARG")]=eij[("ARG","GLY")]=-1.75894
eij[("GLY","LYS")]=eij[("LYS","GLY")]=-1.60783
eij[("GLY","PRO")]=eij[("PRO","GLY")]=-1.42944
eij[("THR","THR")]=-2.27205
eij[("THR","SER")]=eij[("SER","THR")]=-2.44856
eij[("THR","ASN")]=eij[("ASN","THR")]=-3.12074
eij[("THR","GLN")]=eij[("GLN","THR")]=-2.30994
eij[("THR","ASP")]=eij[("ASP","THR")]=-3.37199
eij[("THR","GLU")]=eij[("GLU","THR")]=-2.30207
eij[("THR","HIS")]=eij[("HIS","THR")]=-1.99645
eij[("THR","ARG")]=eij[("ARG","THR")]=-2.12535
eij[("THR","LYS")]=eij[("LYS","THR")]=-2.00369
eij[("THR","PRO")]=eij[("PRO","THR")]=-1.77473
eij[("SER","SER")]=-2.71935
eij[("SER","ASN")]=eij[("ASN","SER")]=-2.88679
eij[("SER","GLN")]=eij[("GLN","SER")]=-2.31760
eij[("SER","ASP")]=eij[("ASP","SER")]=-3.35909
eij[("SER","GLU")]=eij[("GLU","SER")]=-2.12226
eij[("SER","HIS")]=eij[("HIS","SER")]=-2.21046
eij[("SER","ARG")]=eij[("ARG","SER")]=-2.59047
eij[("SER","LYS")]=eij[("LYS","SER")]=-1.96942
eij[("SER","PRO")]=eij[("PRO","SER")]=-1.92283
eij[("ASN","ASN")]=-3.68717
eij[("ASN","GLN")]=eij[("GLN","ASN")]=-2.67566
eij[("ASN","ASP")]=eij[("ASP","ASN")]=-3.64227
eij[("ASN","GLU")]=eij[("GLU","ASN")]=-2.36536
eij[("ASN","HIS")]=eij[("HIS","ASN")]=-2.59329
eij[("ASN","ARG")]=eij[("ARG","ASN")]=-3.00288
eij[("ASN","LYS")]=eij[("LYS","ASN")]=-2.47627
eij[("ASN","PRO")]=eij[("PRO","ASN")]=-2.60250
eij[("GLN","GLN")]=-2.15872
eij[("GLN","ASP")]=eij[("ASP","GLN")]=-1.98926
eij[("GLN","GLU")]=eij[("GLU","GLN")]=-2.27565
eij[("GLN","HIS")]=eij[("HIS","GLN")]=-1.90415
eij[("GLN","ARG")]=eij[("ARG","GLN")]=-2.06515
eij[("GLN","LYS")]=eij[("LYS","GLN")]=-1.60899
eij[("GLN","PRO")]=eij[("PRO","GLN")]=-1.72891
eij[("ASP","ASP")]=-2.14122
eij[("ASP","GLU")]=eij[("GLU","ASP")]=-1.30591
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-2.01424
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-5.52904
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-5.02955
eij[("ASP","PRO")]=eij[("PRO","ASP")]=-2.36134
eij[("GLU","GLU")]=-1.29384
eij[("GLU","HIS")]=eij[("HIS","GLU")]=-1.50644
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-4.42694
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-2.61801
eij[("GLU","PRO")]=eij[("PRO","GLU")]=-1.27329
eij[("HIS","HIS")]=-2.70975
eij[("HIS","ARG")]=eij[("ARG","HIS")]=-3.47470
eij[("HIS","LYS")]=eij[("LYS","HIS")]=-2.32972
eij[("HIS","PRO")]=eij[("PRO","HIS")]=-1.41680
eij[("ARG","ARG")]=-1.83495
eij[("ARG","LYS")]=eij[("LYS","ARG")]=-1.45082
eij[("ARG","PRO")]=eij[("PRO","ARG")]=-1.45422
eij[("LYS","LYS")]=-0.59997
eij[("LYS","PRO")]=eij[("PRO","LYS")]=-1.74852
eij[("PRO","PRO")]=-1.50787


list=['CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR','ALA','GLY','THR','SER','ASN','GLN','ASP','GLU','HIS','ARG','LYS','PRO']

def Compute_ekr(resi):
    ekr=0
    for i in list:
        ekr=ekr+eij[(resi,i)]

    ekr=ekr/20
    return ekr
    

def Compute_err():
    err=0
    for i in list:
        err=err+Compute_ekr(i)

    err=err/20
    return err


ekr={}
ekr['CYS']=Compute_ekr('CYS')
ekr['MET']=Compute_ekr('MET')
ekr['PHE']=Compute_ekr('PHE')
ekr['ILE']=Compute_ekr('ILE')
ekr['LEU']=Compute_ekr('LEU')
ekr['VAL']=Compute_ekr('VAL')
ekr['TRP']=Compute_ekr('TRP')
ekr['TYR']=Compute_ekr('TYR')
ekr['ALA']=Compute_ekr('ALA')
ekr['GLY']=Compute_ekr('GLY')
ekr['THR']=Compute_ekr('THR')
ekr['SER']=Compute_ekr('SER')
ekr['ASN']=Compute_ekr('ASN')
ekr['GLN']=Compute_ekr('GLN')
ekr['ASP']=Compute_ekr('ASP')
ekr['GLU']=Compute_ekr('GLU')
ekr['HIS']=Compute_ekr('HIS')
ekr['ARG']=Compute_ekr('ARG')
ekr['LYS']=Compute_ekr('LYS')
ekr['PRO']=Compute_ekr('PRO')

err=Compute_err()


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
      

def scoreMDW(pdbid, ctype):


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
            with open('./raw_graphv2/'+pdbid+'_CA_output.txt', 'a') as f:
                print(line, file=f)
        if (ctype == 'CB'):
            with open('./raw_graphv2/'+pdbid+'_CB_output.txt', 'a') as f:
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
            #print("resPid:"+resPid+" Etot:"+str(round(Etot,2))+" Ec:"+str(round(Ec,2))+" Ep:"+str(round(Ep,2))+" lines:"+str(num_lines)+"/"+str(total_lines))
            Ec=0
            Ep=0
            rpack=0
            resPid=resAid
    
        Ec=Ec+Eij
        rpack=rpack+1

    print("scoreMDW:",str(round(Etot,2))," ctype:",ctype)
    return round(Etot,2)

