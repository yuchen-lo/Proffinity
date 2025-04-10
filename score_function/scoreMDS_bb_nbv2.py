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
eij[("CYS","CYS")]=-2.17904
eij[("CYS","MET")]=eij[("MET","CYS")]=-2.78557
eij[("CYS","PHE")]=eij[("PHE","CYS")]=-2.07131
eij[("CYS","ILE")]=eij[("ILE","CYS")]=-3.07283
eij[("CYS","LEU")]=eij[("LEU","CYS")]=-2.94957
eij[("CYS","VAL")]=eij[("VAL","CYS")]=-2.91358
eij[("CYS","TRP")]=eij[("TRP","CYS")]=-2.29893
eij[("CYS","TYR")]=eij[("TYR","CYS")]=-2.23920
eij[("CYS","ALA")]=eij[("ALA","CYS")]=-2.44362
eij[("CYS","GLY")]=eij[("GLY","CYS")]=-2.44362
eij[("CYS","THR")]=eij[("THR","CYS")]=-2.25201
eij[("CYS","SER")]=eij[("SER","CYS")]=-1.74890
eij[("CYS","ASN")]=eij[("ASN","CYS")]=-2.14331
eij[("CYS","GLN")]=eij[("GLN","CYS")]=-2.02872
eij[("CYS","ASP")]=eij[("ASP","CYS")]=-0.88802
eij[("CYS","GLU")]=eij[("GLU","CYS")]=-0.93326
eij[("CYS","HIS")]=eij[("HIS","CYS")]=-2.35508
eij[("CYS","ARG")]=eij[("ARG","CYS")]=-2.40396
eij[("CYS","LYS")]=eij[("LYS","CYS")]=-1.69258
eij[("CYS","PRO")]=eij[("PRO","CYS")]=-2.41110
eij[("MET","MET")]=-3.38118
eij[("MET","PHE")]=eij[("PHE","MET")]=-2.94140
eij[("MET","ILE")]=eij[("ILE","MET")]=-3.64020
eij[("MET","LEU")]=eij[("LEU","MET")]=-3.65723
eij[("MET","VAL")]=eij[("VAL","MET")]=-3.61291
eij[("MET","TRP")]=eij[("TRP","MET")]=-3.30376
eij[("MET","TYR")]=eij[("TYR","MET")]=-3.40879
eij[("MET","ALA")]=eij[("ALA","MET")]=-2.46192
eij[("MET","GLY")]=eij[("GLY","MET")]=-2.46192
eij[("MET","THR")]=eij[("THR","MET")]=-2.44840
eij[("MET","SER")]=eij[("SER","MET")]=-2.41442
eij[("MET","ASN")]=eij[("ASN","MET")]=-2.54180
eij[("MET","GLN")]=eij[("GLN","MET")]=-2.50619
eij[("MET","ASP")]=eij[("ASP","MET")]=-0.91324
eij[("MET","GLU")]=eij[("GLU","MET")]=-1.22834
eij[("MET","HIS")]=eij[("HIS","MET")]=-2.64525
eij[("MET","ARG")]=eij[("ARG","MET")]=-3.15663
eij[("MET","LYS")]=eij[("LYS","MET")]=-1.99843
eij[("MET","PRO")]=eij[("PRO","MET")]=-2.85824
eij[("PHE","PHE")]=-3.06682
eij[("PHE","ILE")]=eij[("ILE","PHE")]=-2.96207
eij[("PHE","LEU")]=eij[("LEU","PHE")]=-2.84502
eij[("PHE","VAL")]=eij[("VAL","PHE")]=-2.76432
eij[("PHE","TRP")]=eij[("TRP","PHE")]=-3.56630
eij[("PHE","TYR")]=eij[("TYR","PHE")]=-3.43917
eij[("PHE","ALA")]=eij[("ALA","PHE")]=-2.01230
eij[("PHE","GLY")]=eij[("GLY","PHE")]=-2.01230
eij[("PHE","THR")]=eij[("THR","PHE")]=-1.95279
eij[("PHE","SER")]=eij[("SER","PHE")]=-1.84605
eij[("PHE","ASN")]=eij[("ASN","PHE")]=-2.03255
eij[("PHE","GLN")]=eij[("GLN","PHE")]=-1.79398
eij[("PHE","ASP")]=eij[("ASP","PHE")]=-0.99498
eij[("PHE","GLU")]=eij[("GLU","PHE")]=-1.12033
eij[("PHE","HIS")]=eij[("HIS","PHE")]=-2.37259
eij[("PHE","ARG")]=eij[("ARG","PHE")]=-4.50374
eij[("PHE","LYS")]=eij[("LYS","PHE")]=-2.08061
eij[("PHE","PRO")]=eij[("PRO","PHE")]=-2.29473
eij[("ILE","ILE")]=-4.01318
eij[("ILE","LEU")]=eij[("LEU","ILE")]=-3.76265
eij[("ILE","VAL")]=eij[("VAL","ILE")]=-3.96035
eij[("ILE","TRP")]=eij[("TRP","ILE")]=-2.94878
eij[("ILE","TYR")]=eij[("TYR","ILE")]=-3.37080
eij[("ILE","ALA")]=eij[("ALA","ILE")]=-2.99546
eij[("ILE","GLY")]=eij[("GLY","ILE")]=-2.99553
eij[("ILE","THR")]=eij[("THR","ILE")]=-2.99553
eij[("ILE","SER")]=eij[("SER","ILE")]=-2.76449
eij[("ILE","ASN")]=eij[("ASN","ILE")]=-2.52536
eij[("ILE","GLN")]=eij[("GLN","ILE")]=-3.00896
eij[("ILE","ASP")]=eij[("ASP","ILE")]=-1.01717
eij[("ILE","GLU")]=eij[("GLU","ILE")]=-1.40244
eij[("ILE","HIS")]=eij[("HIS","ILE")]=-2.97079
eij[("ILE","ARG")]=eij[("ARG","ILE")]=-3.89904
eij[("ILE","LYS")]=eij[("LYS","ILE")]=-2.13992
eij[("ILE","PRO")]=eij[("PRO","ILE")]=-3.42980
eij[("LEU","LEU")]=-3.82304
eij[("LEU","VAL")]=eij[("VAL","LEU")]=-3.86267
eij[("LEU","TRP")]=eij[("TRP","LEU")]=-3.02414
eij[("LEU","TYR")]=eij[("TYR","LEU")]=-3.15649
eij[("LEU","ALA")]=eij[("ALA","LEU")]=-2.83483
eij[("LEU","GLY")]=eij[("GLY","LEU")]=-2.83483
eij[("LEU","THR")]=eij[("THR","LEU")]=-2.91129
eij[("LEU","SER")]=eij[("SER","LEU")]=-2.48819
eij[("LEU","ASN")]=eij[("ASN","LEU")]=-2.77837
eij[("LEU","GLN")]=eij[("GLN","LEU")]=-3.07365
eij[("LEU","ASP")]=eij[("ASP","LEU")]=-0.92955
eij[("LEU","GLU")]=eij[("GLU","LEU")]=-1.67601
eij[("LEU","HIS")]=eij[("HIS","LEU")]=-3.17764
eij[("LEU","ARG")]=eij[("ARG","LEU")]=-3.66126
eij[("LEU","LYS")]=eij[("LYS","LEU")]=-2.07002
eij[("LEU","PRO")]=eij[("PRO","LEU")]=-3.35522
eij[("VAL","VAL")]=-3.79238
eij[("VAL","TRP")]=eij[("TRP","VAL")]=-2.73813
eij[("VAL","TYR")]=eij[("TYR","VAL")]=-2.96476
eij[("VAL","ALA")]=eij[("ALA","VAL")]=-3.11309
eij[("VAL","GLY")]=eij[("GLY","VAL")]=-3.11309
eij[("VAL","THR")]=eij[("THR","VAL")]=-2.86107
eij[("VAL","SER")]=eij[("SER","VAL")]=-2.86595
eij[("VAL","ASN")]=eij[("ASN","VAL")]=-2.72058
eij[("VAL","GLN")]=eij[("GLN","VAL")]=-2.88787
eij[("VAL","ASP")]=eij[("ASP","VAL")]=-1.03742
eij[("VAL","GLU")]=eij[("GLU","VAL")]=-1.33513
eij[("VAL","HIS")]=eij[("HIS","VAL")]=-2.62852
eij[("VAL","ARG")]=eij[("ARG","VAL")]=-3.13778
eij[("VAL","LYS")]=eij[("LYS","VAL")]=-2.37129
eij[("VAL","PRO")]=eij[("PRO","VAL")]=-3.15348
eij[("TRP","TRP")]=-4.39703
eij[("TRP","TYR")]=eij[("TYR","TRP")]=-3.58668
eij[("TRP","ALA")]=eij[("ALA","TRP")]=-1.47902
eij[("TRP","GLY")]=eij[("GLY","TRP")]=-1.47902
eij[("TRP","THR")]=eij[("THR","TRP")]=-2.10644
eij[("TRP","SER")]=eij[("SER","TRP")]=-1.55757
eij[("TRP","ASN")]=eij[("ASN","TRP")]=-2.34506
eij[("TRP","GLN")]=eij[("GLN","TRP")]=-3.05325
eij[("TRP","ASP")]=eij[("ASP","TRP")]=-1.02100
eij[("TRP","GLU")]=eij[("GLU","TRP")]=-1.02386
eij[("TRP","HIS")]=eij[("HIS","TRP")]=-2.93602
eij[("TRP","ARG")]=eij[("ARG","TRP")]=-7.03968
eij[("TRP","LYS")]=eij[("LYS","TRP")]=-4.77708
eij[("TRP","PRO")]=eij[("PRO","TRP")]=-3.08244
eij[("TYR","TYR")]=-3.63958
eij[("TYR","ALA")]=eij[("ALA","TYR")]=-1.74571
eij[("TYR","GLY")]=eij[("GLY","TYR")]=-1.74571
eij[("TYR","THR")]=eij[("THR","TYR")]=-2.17124
eij[("TYR","SER")]=eij[("SER","TYR")]=-1.82703
eij[("TYR","ASN")]=eij[("ASN","TYR")]=-2.51450
eij[("TYR","GLN")]=eij[("GLN","TYR")]=-2.80451
eij[("TYR","ASP")]=eij[("ASP","TYR")]=-1.05010
eij[("TYR","GLU")]=eij[("GLU","TYR")]=-1.26683
eij[("TYR","HIS")]=eij[("HIS","TYR")]=-3.04387
eij[("TYR","ARG")]=eij[("ARG","TYR")]=-5.72831
eij[("TYR","LYS")]=eij[("LYS","TYR")]=-3.47643
eij[("TYR","PRO")]=eij[("PRO","TYR")]=-2.99753
eij[("ALA","ALA")]=-3.01226
eij[("ALA","GLY")]=eij[("GLY","ALA")]=-3.01226
eij[("ALA","THR")]=eij[("THR","ALA")]=-2.39220
eij[("ALA","SER")]=eij[("SER","ALA")]=-2.09699
eij[("ALA","ASN")]=eij[("ASN","ALA")]=-2.06177
eij[("ALA","GLN")]=eij[("GLN","ALA")]=-2.00413
eij[("ALA","ASP")]=eij[("ASP","ALA")]=-0.50172
eij[("ALA","GLU")]=eij[("GLU","ALA")]=-0.84664
eij[("ALA","HIS")]=eij[("HIS","ALA")]=-1.94311
eij[("ALA","ARG")]=eij[("ARG","ALA")]=-1.43226
eij[("ALA","LYS")]=eij[("LYS","ALA")]=-1.70043
eij[("ALA","PRO")]=eij[("PRO","ALA")]=-2.80274
eij[("GLY","GLY")]=-3.01226
eij[("GLY","THR")]=eij[("THR","GLY")]=-2.39220
eij[("GLY","SER")]=eij[("SER","GLY")]=-2.09699
eij[("GLY","ASN")]=eij[("ASN","GLY")]=-2.06177
eij[("GLY","GLN")]=eij[("GLN","GLY")]=-2.00413
eij[("GLY","ASP")]=eij[("ASP","GLY")]=-0.50172
eij[("GLY","GLU")]=eij[("GLU","GLY")]=-0.84664
eij[("GLY","HIS")]=eij[("HIS","GLY")]=-1.94311
eij[("GLY","ARG")]=eij[("ARG","GLY")]=-1.43226
eij[("GLY","LYS")]=eij[("LYS","GLY")]=-1.70043
eij[("GLY","PRO")]=eij[("PRO","GLY")]=-2.80274
eij[("THR","THR")]=-2.17639
eij[("THR","SER")]=eij[("SER","THR")]=-2.15427
eij[("THR","ASN")]=eij[("ASN","THR")]=-2.02294
eij[("THR","GLN")]=eij[("GLN","THR")]=-2.33920
eij[("THR","ASP")]=eij[("ASP","THR")]=-0.63305
eij[("THR","GLU")]=eij[("GLU","THR")]=-1.13169
eij[("THR","HIS")]=eij[("HIS","THR")]=-1.99020
eij[("THR","ARG")]=eij[("ARG","THR")]=-1.89405
eij[("THR","LYS")]=eij[("LYS","THR")]=-1.59374
eij[("THR","PRO")]=eij[("PRO","THR")]=-2.45346
eij[("SER","SER")]=-1.80919
eij[("SER","ASN")]=eij[("ASN","SER")]=-1.97778
eij[("SER","GLN")]=eij[("GLN","SER")]=-2.37238
eij[("SER","ASP")]=eij[("ASP","SER")]=-0.57952
eij[("SER","GLU")]=eij[("GLU","SER")]=-1.32432
eij[("SER","HIS")]=eij[("HIS","SER")]=-2.04149
eij[("SER","ARG")]=eij[("ARG","SER")]=-2.07649
eij[("SER","LYS")]=eij[("LYS","SER")]=-1.48734
eij[("SER","PRO")]=eij[("PRO","SER")]=-2.14859
eij[("ASN","ASN")]=-2.19598
eij[("ASN","GLN")]=eij[("GLN","ASN")]=-2.28060
eij[("ASN","ASP")]=eij[("ASP","ASN")]=-1.04591
eij[("ASN","GLU")]=eij[("GLU","ASN")]=-1.29822
eij[("ASN","HIS")]=eij[("HIS","ASN")]=-1.87940
eij[("ASN","ARG")]=eij[("ARG","ASN")]=-1.92194
eij[("ASN","LYS")]=eij[("LYS","ASN")]=-1.33357
eij[("ASN","PRO")]=eij[("PRO","ASN")]=-2.25369
eij[("GLN","GLN")]=-2.23722
eij[("GLN","ASP")]=eij[("ASP","GLN")]=-1.29546
eij[("GLN","GLU")]=eij[("GLU","GLN")]=-1.09240
eij[("GLN","HIS")]=eij[("HIS","GLN")]=-1.96400
eij[("GLN","ARG")]=eij[("ARG","GLN")]=-2.26601
eij[("GLN","LYS")]=eij[("LYS","GLN")]=-1.62957
eij[("GLN","PRO")]=eij[("PRO","GLN")]=-2.79666
eij[("ASP","ASP")]=-0.61225
eij[("ASP","GLU")]=eij[("GLU","ASP")]=-0.65366
eij[("ASP","HIS")]=eij[("HIS","ASP")]=-1.41218
eij[("ASP","ARG")]=eij[("ARG","ASP")]=-7.79878
eij[("ASP","LYS")]=eij[("LYS","ASP")]=-3.25082
eij[("ASP","PRO")]=eij[("PRO","ASP")]=-1.06597
eij[("GLU","GLU")]=-0.80181
eij[("GLU","HIS")]=eij[("HIS","GLU")]=-1.61068
eij[("GLU","ARG")]=eij[("ARG","GLU")]=-7.45705
eij[("GLU","LYS")]=eij[("LYS","GLU")]=-2.76222
eij[("GLU","PRO")]=eij[("PRO","GLU")]=-1.63794
eij[("HIS","HIS")]=-2.83985
eij[("HIS","ARG")]=eij[("ARG","HIS")]=-3.41747
eij[("HIS","LYS")]=eij[("LYS","HIS")]=-1.94223
eij[("HIS","PRO")]=eij[("PRO","HIS")]=-2.57331
eij[("ARG","ARG")]=-1.39272
eij[("ARG","LYS")]=eij[("LYS","ARG")]=-0.90115
eij[("ARG","PRO")]=eij[("PRO","ARG")]=-2.74562
eij[("LYS","LYS")]=-0.39426
eij[("LYS","PRO")]=eij[("PRO","LYS")]=-1.88178
eij[("PRO","PRO")]=-3.02585


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
      

def scoreMDS_BB_NB(pdbid, ctype):

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
            #print("resPid:"+resPid+" Etot:"+str(round(Etot/100,2))+" Ec:"+str(round(Ec,2))+" lines:"+str(num_lines)+"/"+str(total_lines))
            Ec=0
            Ep=0
            rpack=0
            resPid=resAid
            
        Ec=Ec+Eij
        rpack=rpack+1

    print("scoreMDS_BB_NB:",str(round(Etot/100,2))," ctype:",ctype)
    return round(Etot/100,2)
