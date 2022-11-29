from general_utils.mrc_uilts import get_mrc_level, get_mrc_to_pdb_aux
from general_utils.pdb_utils import align_pdb_file_1_in_2, align_tmaling

mrc_path = "/home/lcastillo98/Downloads/Datos_de_Prueba/1yfq_A.mrc"
pdb_path = "/home/lcastillo98/Downloads/Datos_de_Prueba/1yfq_CA.pdb"


simulate_pdb_path = "/home/lcastillo98/Downloads/Datos_de_Prueba/1yfq_CA_simulate.pdb"
simulate_good_pdb_path = "/home/lcastillo98/Downloads/Datos_de_Prueba/1yfq_CA_simulate_good.pdb"

import pymol2

#get_mrc_to_pdb_aux(9.94, mrc_path, simulate_pdb_path)

result = align_tmaling(simulate_pdb_path, pdb_path)
print("TM_score_normalized_Chain_1:", result.TM_score_normalized_Chain_1)
print("TM_score_normalized_Chain_2:", result.TM_score_normalized_Chain_2)
print("RMS:", result.RMSD)
print("aligned_length:", result.aligned_length)
print("Seq_ID:", result.Seq_ID)
print("length_of_Chain_1:", result.length_of_Chain_1)
print("length_of_Chain_2:", result.length_of_Chain_2)

print("\n\n")
result = align_tmaling(simulate_good_pdb_path, pdb_path)
print("TM_score_normalized_Chain_1:", result.TM_score_normalized_Chain_1)
print("TM_score_normalized_Chain_2:", result.TM_score_normalized_Chain_2)
print("RMS:", result.RMSD)
print("aligned_length:", result.aligned_length)
print("Seq_ID:", result.Seq_ID)
print("length_of_Chain_1:", result.length_of_Chain_1)
print("length_of_Chain_2:", result.length_of_Chain_2)



# from tmtools.io import get_structure, get_residue_data
# from tmtools import tm_align
#
# s1 = get_structure(pdb_path)
# s2 = get_structure(simulate_pdb_path)
#
# coords1, seq1 = get_residue_data(s1)
# coords2, seq2 = get_residue_data(s2)
#
# res = tm_align(coords1, coords2, seq1, seq2)
# print(res.tm_norm_chain1)
# print(res.tm_norm_chain2)
# print(res.t)
# print(res.u)


#level = get_mrc_level(mrc_path, pdb_path)
#print(level)

#aling_result = align_pdb_file_1_in_2(simulate_pdb_path, pdb_path)
#print(aling_result.RMSDAfterRefinement)

