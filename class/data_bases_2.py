from Bio import Entrez, SeqIO
import pickle

# Poner correo electrónico para entrar a Entrez
Entrez.email = "dgoretti@lcg.unam.mx"

# EGQUERY
termino = "(Aedes[Title] OR Aedes[All Fields])AND((RNA-Seq[Title] OR transcriptomic[Title]) OR (transcriptome[Title] " \
          "OR sequencing[Title]))"
handle = Entrez.egquery(term=termino)
record = Entrez.read(handle)

# ESPELL
handle = Entrez.espell(term="biopythooon")
record1 = Entrez.read(handle)

# EFETCH
# db = "nuccore" también es válido
ids = "1919569438, 1919569357, 1251949171"
handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
# usamos parse porque hay más de un read
record2 = SeqIO.parse(handle, "genbank")

# ESUMMARY
handle = Entrez.esummary(db="taxonomy", id="9913,30521")
record3 = Entrez.read(handle)
print(len(record3))

# tamaño del record deseado - tamaño esummary
print(len(pickle.dumps(record)))

# EFETCH Y ARCHIVOS TEXTO
out_handle = open("files/prueba.fasta", "w")
fetch_handle = Entrez.efetch(db="nucleotide", id="1919569438, 1919569357, 1251949171", rettype="fasta", retmode="text")
# Usar handle.read()
data = fetch_handle.read()
# Cerrar handle
fetch_handle.close()
# Escribir archivo
out_handle.write(data)
# Cerrar archivo
out_handle.close()


# EJERCICIO 2 ---------------------------------------------------------------------------

# Notoryctes typhlops----------------------------------------------------
# PRIMERA PARTE: esearch para buscar 1er organismo en taxonomy
handle = Entrez.esearch(db="Taxonomy", term="Notoryctes typhlops")
record4 = Entrez.read(handle)
# obtenemos ID de taxonomía
print(record4["IdList"])
# Guardar ID
id_taxo = record4["IdList"]
# SEGUNDA PARTE: efetch para obtener archivo de taxonomía
handle = Entrez.efetch(db="Taxonomy", id=id_taxo, retmode="xml")
notoryctes = Entrez.read(handle)
# Checar información
print(notoryctes[0].keys())
print(handle.url)

# Ver linaje
print(notoryctes[0]["Lineage"])

# Chrysochloris asiatica -----------------------------------------------
# PRIMERA PARTE
handle = Entrez.esearch(db="Taxonomy", term="Chrysochloris asiatica")
record = Entrez.read(handle)
id_taxo = record["IdList"][0]
# SEGUNDA PARTE
handle = Entrez.efetch(db="Taxonomy", id=id_taxo, retmode="xml")
chryso = Entrez.read(handle)
print(chryso[0]["Lineage"])
print(handle.url)

# Ver linaje
print(chryso[0]["Lineage"])


# Checar diferencia entre linajes de especies
topo1 = notoryctes[0]["Lineage"].split(";")
topo2 = chryso[0]["Lineage"].split(";")


def comparar(org1, org2):
    i = -1
    for linaje1, linaje2 in zip(org1, org2):
        i += 1
        if linaje1 != linaje2:
            print(i)
            print("Organismo 1:", linaje1, ";", "Organismo 2:", linaje2)
            break


print(comparar(topo1, topo2))


# ELINK -> ver en einfo para ver su relación
ids = "15718680"
# Buscar con elink los ids de protein en la base de datos de gene
record = Entrez.read(Entrez.elink(dbfrom="protein", id=ids, db='gene'))
# Visualizar record
print(record[0])

# Obtener citas
pmid = "32703847"  # pubmed id
results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid))

# Guardar links de PMC
pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
print(pmc_ids)

# IDs para pubmed:
results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_refs_pubmed", id=",".join(pmc_ids)))
# Guardar links
pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
print(pubmed_ids)
