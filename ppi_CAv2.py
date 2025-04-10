from Bio import PDB
from Bio.PDB import Selection, NeighborSearch
from Bio.PDB import PICIO, PDBIO
from typing import TypedDict, Dict, Tuple
from math import sqrt
from Bio.PDB import PDBParser, Select
import Bio
import pandas as pd

class CloseSelect_7(Select):
    def accept_atom(self, atom):
        if atom.get_serial_number() in close_atoms_number_7:
            return 1
        else:
            return 0

def ppi_CA(pdbid, cutoff):


    print(pdbid+" is the target protein")
    print(str(cutoff)+" is the cutoff")

    p = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    structure = p.get_structure("X", "./model/"+str(pdbid)+".pdb")

    list_all=[]
    list_structure=[]
    
    list_ca_bb_b_all=[]
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    s=structure[0]
                    if (atom.get_id()=='CA' and residue.get_id()[0] == " "):
                        target_atom = s[chain.get_id()][residue.get_id()]['CA']
                    
                    
                        list=[chain.get_id(), residue.get_resname(),residue.get_id()[1],target_atom.coord]
                        list_all.append(list)
                        list_structure.append(list)
                    

        for i in list_structure:
            for j in list_structure:
                if (i[0] == j[0]):
                    if (i[2]==j[2]+1 or i[2]==j[2]-1):
                        dist=sqrt((i[3][0]-j[3][0])**2+(i[3][1]-j[3][1])**2+(i[3][2]-j[3][2])**2)
                        list=[i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2]]
                        list_ca_bb_b_all.append(list)
                    #with open(args.pdbid+"_CA_output_bb_b.txt", "a") as f:
                        #print(i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2], file=f)  


    list_ca_all=[]
    list_ca_bb_nb_all=[]
    
    for i in list_all:
       for j in list_all:
          if (i[0] != j[0]):
              #print("not eq")
              dist=sqrt((i[3][0]-j[3][0])**2+(i[3][1]-j[3][1])**2+(i[3][2]-j[3][2])**2)
              if (dist<=float(cutoff)):
                  list=[i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2]]
                  list_ca_all.append(list)
                	   #with open(args.pdbid+"_CA_output.txt", "a") as f:
                	      #print(i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2], file=f)

          elif (i[0] == j[0]):
              #print("eq")
              dist=sqrt((i[3][0]-j[3][0])**2+(i[3][1]-j[3][1])**2+(i[3][2]-j[3][2])**2)
              if ((i[2]>j[2]+1 or i[2]<j[2]-1) and dist<=float(cutoff)):
                  list=[i[0],i[1],i[2],j[0],j[1],j[2],dist, i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2]]
                  list_ca_bb_nb_all.append(list)
                 
               #with open(args.pdbid+"_CA_output_bb_nb.txt", "a") as f:
                   #print(i[0],i[1],i[2],j[0],j[1],j[2],dist, i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2],file=f)


    df_ca_all=pd.DataFrame(list_ca_all)
    #df_ca_all=df_ca_all.drop_duplicates()
    df_ca_all.to_csv("./raw_graphv2/"+str(pdbid)+"_CA_output.txt", index=False, header=False)
    df_ca_bb_b_all=pd.DataFrame(list_ca_bb_b_all)
    df_ca_bb_b_all=df_ca_bb_b_all.drop_duplicates()
    df_ca_bb_b_all.to_csv("./raw_graphv2/"+str(pdbid)+"_CA_output_bb_b.txt", index=False, header=False)
    df_ca_bb_nb_all=pd.DataFrame(list_ca_bb_nb_all)
    df_ca_bb_nb_all=df_ca_bb_nb_all.drop_duplicates()
    df_ca_bb_nb_all.to_csv("./raw_graphv2/"+str(pdbid)+"_CA_output_bb_nb.txt", index=False, header=False)

    #return df_ca_all, df_ca_bb_nb_all, df_ca_bb_b_all
    return


