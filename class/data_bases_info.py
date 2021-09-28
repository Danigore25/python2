from Bio import Entrez
from pprint import pprint

# Correo
Entrez.email = "dgoretti@lcg.unam.mx"

# handle con einfo en Genome
handle = Entrez.einfo(db="genome")
record = Entrez.read(handle)
print(record["DbInfo"]["Description"])
print(record["DbInfo"]["FieldList"])
print(handle.url)
handle.close()

# EJERCICIO 1

i = -1
for field in record["DbInfo"]["FieldList"]:
    i += 1
    if field["Name"] == "ORGN":
        print(i, field["Name"], field["Description"])

# Acceder a índice de Organism obtenido del anterior paso
print(record["DbInfo"]["FieldList"][3]["Description"])

# handle con esearch
handle = Entrez.esearch(db="pubmed", term="biopython")
record2 = Entrez.read(handle)
print(record2["Count"])
handle.close()

# Ejemplo buscando autora
handle = Entrez.esearch(db="pubmed", term="Valeria Mateo-Estrada", field="AUTH")
record3 = Entrez.read(handle)
print(record3["IdList"])
handle.close()
