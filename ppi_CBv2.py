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

# Define the parser
def ppi_CB(pdbid, cutoff):

    print(pdbid+" is the target protein")
    print(str(cutoff) +" is the cutoff")

    p = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    structure = p.get_structure("X", "./model/"+str(pdbid)+".pdb")

    list_all=[]

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    s=structure[0]
                    if (atom.get_id()=='CB' and residue.get_id()[0] == " "):
                        target_atom = s[chain.get_id()][residue.get_id()]['CB']
                        list=[chain.get_id(),residue.get_resname(),residue.get_id()[1],target_atom.coord]
                        list_all.append(list)
                    elif (residue.get_resname()=='GLY' and atom.get_id()=='CA'):
                        target_atom = s[chain.get_id()][residue.get_id()]['CA']
                        list=[chain.get_id(),residue.get_resname(),residue.get_id()[1],target_atom.coord]
                        list_all.append(list)
    list_cb_all=[]
    list_cb_bb_all=[]
    
    for i in list_all:
       for j in list_all:
          if (i[0] != j[0]):
            	dist=sqrt((i[3][0]-j[3][0])**2+(i[3][1]-j[3][1])**2+(i[3][2]-j[3][2])**2)
            	if (dist<=float(cutoff)):
                    list=[i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2]]
                    list_cb_all.append(list)
                	   #with open(pdbid+"_CB_output.txt", "a") as f:
                	      #print(i[0],i[1],i[2],j[0],j[1],j[2],dist, i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2], file=f)

          elif (i[0] == j[0]):
              dist=sqrt((i[3][0]-j[3][0])**2+(i[3][1]-j[3][1])**2+(i[3][2]-j[3][2])**2)
              if (dist>0 and dist<=float(cutoff)):
                  list=[i[0],i[1],i[2],j[0],j[1],j[2],dist,i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2]]
                  list_cb_bb_all.append(list)
                #with open(pdbid+"_CB_output_bb.txt", "a") as f:
                   #print(i[0],i[1],i[2],j[0],j[1],j[2],dist, i[3][0],i[3][1],i[3][2],j[3][0],j[3][1],j[3][2], file=f)
    
    df_cb_all=pd.DataFrame(list_cb_all)
    #df_cb_all=df_cb_all.drop_duplicates()
    df_cb_all.to_csv("./raw_graphv2/"+str(pdbid)+"_CB_output.txt", index=False, header=False)
    df_cb_bb_all=pd.DataFrame(list_cb_bb_all)
    df_cb_bb_all=df_cb_bb_all.drop_duplicates()
    df_cb_bb_all.to_csv("./raw_graphv2/"+str(pdbid)+"_CB_output_bb.txt", index=False, header=False)
    
    #return df_cb_all, df_cb_bb_all
    return


