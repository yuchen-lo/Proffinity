from modeller import *

env = Environ()
env.io.atom_files_directory = ['./PDBs/']
aln = Alignment(env)
mdl = Model(env, file='1CSE', model_segment=('FIRST:E','LAST:I'))
aln.append_model(mdl, align_codes='1CSE-E-I', atom_files='1CSE.pdb')
aln.append(file='seq1_out.ali', align_codes='1CSE')
aln.align2d(max_gap_length=50)
aln.write(file='template-target.ali', alignment_format='PIR')
aln.write(file='template-target.pap', alignment_format='PAP')
