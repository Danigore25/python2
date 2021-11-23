import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from Bio.Seq import Seq
from Bio import SeqIO
import csv
import re


count_matrix = np.array([[3, 3, 0], [0, 0, 1], [1, 1, 0], [0, 0, 1], [1, 0, 4], [1, 2, 0]])

count_matrix = np.where(count_matrix > 0, 1, 0)
print(count_matrix)

# Producto punto de matrices
# print("producto punto gen1vsgen1:", np.dot(count_matrix.T[0], count_matrix.T[0]))
#
# #### Gen 1 vs gen 3
# print(np.vstack((count_matrix.T[0],count_matrix.T[2])) ) # para visualizar
#
# print("producto punto gen1 vs gen3:", np.dot(count_matrix.T[0], count_matrix.T[2]))

# Hacer multiplicación de matrices
expresion = np.matmul(count_matrix, count_matrix)
print('Expresión: ')
print(expresion)

# Plotear
plt.imshow(expresion)
plt.colorbar()
plt.show()

# Grafo
G = nx.DiGraph(expresion)
nx.draw(G, node_size=900,  with_labels=True)
plt.show()

# EJERCICIO 1


def secuencia_aleatoria(tamano=100, p=[0.5, 0.5, 0.5, 0.5], seed=None):
    np.random.seed(seed)  # posibilidad de reproducibilidad
    dna = list("ATGC")
    # secuencia random con distribucion p
    secuencia2 = Seq(''.join(np.random.choice(dna, tamano, p)))
    return secuencia2
# probemos


secuencia3 = secuencia_aleatoria(25, p=[0.1, 0.2, 0.4, 0.3])
print(secuencia3)

# Hacer secuencias aleatorias sin numpy


def seq_aleatorias(tamano=100, p=[0.5, 0.5, 0.5, 0.5], seed=None,
                   num_seq=10, archivo="files/SeqAleatorias.fasta"):
    # escribamos archivo de secuencias!!
    with open(archivo, 'w') as file:
        for i in range(num_seq):
            sec = secuencia_aleatoria(tamano, p, seed)
            file.write(">Seq" + str(i) + "\n")
            file.write(str(sec))
            file.write("\n")
    file.close()

# Con fasta y
# from Bio import SeqIO
# import csv
# import re


# funcion para busqueda de patrones en archivo fasta


def busqueda_patron(fasta="files/SeqAleatorias.fasta", output="files/ejercicioNP.csv"):
    tf_interes_mod = {"TF_1": 'ATG[GG|TAG]', "TF_2": 'T[TC|AA]GAAT',
                      "TF_3": "GTATGCGGGG", "TF_4": "TAT[GT]CC",
                      "TF_5": "TATATA[GT|TG]"}
    outfile = open(output, "w")
    writer = csv.writer(outfile)
    for rec in SeqIO.parse(fasta, "fasta"):
        # secuencia por analizar
        seq = rec.seq
        # array para guardar numero de coindicencias
        tf_counts = np.empty(0, dtype="int")
        for tf in tf_interes_mod.values():
            counts = len(re.findall(tf, str(seq)))
            tf_counts = np.append(tf_counts, [counts])
        # escribir counts encontrados (iterable)
        writer.writerows([tf_counts])
    outfile.close()


###############################################################
# generar fasta con 10 secuencias distribucion p

print(seq_aleatorias(p=[0.3, 0.2, 0.3, 0.2]))

# busqueda de patrones
print(busqueda_patron())

# matriz de cuentas
count_matrix = np.loadtxt("files/ejercicioNP.csv", delimiter=",", dtype="int")
print(count_matrix)


# suma de los ejes 0
print(np.sum(count_matrix, axis=0))
# suma del eje 1
print(np.sum(count_matrix, axis=1))

# matriz binaria (1 y 0)
matrix_uno = np.where(count_matrix > 0, 1, 0)

# multiplicacion de matrices
coexpresion = np.matmul(matrix_uno.T, matrix_uno)
print(coexpresion)

# plot de matriz de coexistencia de TFs
plt.imshow(coexpresion)
plt.colorbar()
plt.show()

# Grafo de matriz coexpresion
G = nx.DiGraph(coexpresion)
nx.draw(G, node_size=900,  with_labels=True)
plt.show()
