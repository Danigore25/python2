# CLASE 3 DE BIOPYTHON

# EJERCICIO 1

from Bio import PDB
parser = PDB.PDBParser(QUIET=True)
struc = parser.get_structure("prot_1fat", "../docs/1fat.pdb")
print("Llaves: ", struc.header.keys())
print("Método de estructura: ", struc.header['structure_method'])
print("Resolución: ", struc.header['resolution'])

model = struc[0]
print(model.child_dict)

chain1 = model['A']

for chain in model:
    print(chain)
    print(chain.child_list)

# for residue in chain1:
#    print(residue)

residue1 = chain[1]
print("Nombre del primer residuo: ", residue1.get_resname())
print("ID del primer residuo: ", residue1.get_id()[1])

# ARCHIVO 1kcw.pdb

# obj_struc = parser.get_structure("prot_1kcw", "../docs/1kcw.pdb")
#
# # Guardar las cisteínas
#
# cisteinas = []
#
# for modelo in obj_struc:
#     for cadena in modelo:
#         if cadena.id == 'A':
#             for residuo in cadena:
#                 if residuo.get_resname() == 'CYS':
#                     cisteinas.append(residuo)
#
#
# print(len(cisteinas))
