from modeller import *

env = Environ()
env.io.atom_files_directory = ['./PDBs/']
aln = Alignment(env)
mdl = Model(env, file='pdb_id', model_segment=('FIRST:chain_1','LAST:chain_2'))
aln.append_model(mdl, align_codes='pdb_tag', atom_files='pdbfile')
aln.append(file='seqX_out.ali', align_codes='pdb_id')
aln.align2d(max_gap_length=50)
aln.write(file='count.ali', alignment_format='PIR')
aln.write(file='count.pap', alignment_format='PAP')
