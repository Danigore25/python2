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

residue2 = chain[22]

for atom in residue2:
    print(atom)

atom_ca = residue2['CB']

print(atom_ca)

print(atom_ca.coord)

print(atom_ca.element)
print(atom_ca.id)


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
# EJERCICIO 3
#
# for atom2 in cisteinas[0]:
#     print(atom2.element, atom2.id)
#
# EJERCICIO 4
#
# pares = []
# for cisteina_1 in cisteinas:
#     for cisteina_2 in cisteinas:
#         if not (cisteina_1 == cisteina_2):
#             if (cisteina_1['SG'] - cisteina_2['SG']) < 8:
#                 print(cisteina_1.coord, cisteina_2.coord)
#                 print(cisteina_1.id, cisteina2.id)
#
# print(len(cisteinas))
