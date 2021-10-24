from Bio import Entrez, SeqIO
from Bio import ExPASy, SwissProt
from Bio import SeqRecord, Seq
from Bio.ExPASy import Prosite, Prodoc, ScanProsite

Entrez.email = "dgoretti@lcg.unam.mx"

# Buscar organismo y gen en la base de datos de proteínas
handle = Entrez.esearch(db="protein", term="Aedes aegypti[Orgn] AND APY[Gene]")
especie = Entrez.read(handle)
handle.close()
print(especie['IdList'])

# Extraer desde la base de datos de protein
handle = Entrez.efetch(db="protein", id="193806340", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
print(record.annotations['db_source'])
db_source = record.annotations['db_source']
apy_prot = record.annotations['accessions']

# Ver si hay algún gen en UniProt - Hacer parseo

# Buscar información
handle = Entrez.esearch(db="gene", term="Aedes aegypti[Orgn] AND apy[Gene]")
record = Entrez.read(handle)
handle = Entrez.efetch(db="gene", id=record['IdList'][0], retmode="xml")
record1 = Entrez.read(handle, "genbank")
print(record1[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn'])

# Buscar en uniprot
sinonimos = record1[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
print([uniprot for uniprot in sinonimos if '_' in uniprot])

# EJERCICIO 1
# Buscar el accession de la base de datos de UniProtKB del gen DEFA del mosquito Aedes aegypti

handle = Entrez.esearch(db="protein", term="Aedes aegypti[Orgn] AND DEFA[Gene]")
record = Entrez.read(handle)
print(record["IdList"])

handle = Entrez.efetch(db="protein", id=record["IdList"][0], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
print(record.annotations["db_source"])
DEFA_prot = record.annotations['accessions']
db_source = record.annotations["db_source"]
print(db_source)


# EXPASY -----------------------------------------------------------------------------------------------

print(apy_prot[0])
handle = ExPASy.get_sprot_raw(apy_prot[0])
print(handle.url)

record = SwissProt.read(handle)
# Imprimir llaves (data_class, molecule_type, accessions, created, sequence_update, etc)
print(record.__dict__.keys())

print(db_source.split(';')[0])
# 'UniProtKB: locus APY_AEDAE, accession P50635'
handle = ExPASy.get_sprot_raw('APY_AEDAE')
record = SwissProt.read(handle)
handle.close()
print(record.entry_name)
print(record.data_class)
print(record.organism)
print(record.sequence[:10])

# EJERCICIO 2
# Con el accession de la base de datos de UniProtKB del gen DEFA del mosquito Aedes aegypti , obtener el archivo de
# UniProt e imprimir: Fecha de creación, Cuándo actualizaron su anotación, ID de taxonomía

handle = ExPASy.get_sprot_raw('P91793')
record = SwissProt.read(handle)
print(handle.url)

# Fecha de creacion, actualización de la anotación, id de taxonomía
print(record.__dict__.keys())
print(record.cross_references)
print(record.created)
print(record.annotation_update)
print(record.taxonomy_id)
print(record.comments[1])
print(record.references[0])


for reference in record.cross_references:
    if 'PROSITE' in reference:
        print(reference[1])


# Desde un archivo descargado con un solo record -----------------------------------

# handle = open("./files/clase_4/O23729.txt")
# record = SwissProt.read(handle)
# print(record.gene_name)
# print(record.sequence)
# print(record.__dict__.keys())

# MULTIPLES RECORD -----------------------------------------------------------------------

# handle = open("./files/clase_4/mini_uniprot.dat")
# descriptions = [record.description for record in SwissProt.parse(handle)]
# # for record in SwissProt.parse(handle):
# #     descriptions.append(record.description)
# print(descriptions[:5])

# CREAR OBJETO SEQRECORD -------------------------------------------------------------------

# help(Bio.SeqRecord.SeqRecord)

# seqRec = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(record.sequence),
#                                  id=record.entry_name,
#                                  name=record.organism,
#                                  description=record.description)
# print(seqRec.__dict__)
# print(seqRec.format('fasta'))

# from Bio.SeqUtils import seq3, seq1, molecular_weight
# prot = seqRec.seq
# # Convertir a código de 3 letras
# seq3(prot)[0:12]

# # Convertir a código de 1 letra
# seq1( seq3(prot)[0:12] )

# molecular_weight(prot, seq_type='protein')

# EJERCICIO 3
# Con el archivo de UniProtKB del gen DEFA del mosquito Aedes aegypti, crear un objeto SeqRecord

objeto_SeqRecord = SeqRecord.SeqRecord(seq=Seq.Seq(record.sequence), id=record.entry_name, name=record.organism,
                                       description=record.description)

print(objeto_SeqRecord.__dict__)

# # Se puede leer con SeqIO pero hay pérdida de información ya que creamos un objeto SeqRecord
# handle = ExPASy.get_sprot_raw('P91793')

# record = SeqIO.read(handle, 'swiss')

# print(record.__dict__.keys())

# PROSITE --------------------------------------------------------------------------------------

handle = ExPASy.get_prosite_raw("PS00785")
record = Prosite.read(handle)
print(record.name)
print(record.type)
print(record.pattern)
print(record.rules)

handle = ExPASy.get_prosite_raw(record.pdoc)
record = Prodoc.read(handle)
print(record.text)

# SCAN PROSITE

sequence = "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFTCRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLS" \
           "ECKTGNGKNYRGTMSKTKN"
handle = ScanProsite.scan(seq=sequence)
result = ScanProsite.read(handle)
print(type(result))
print(result[0])

handle = ExPASy.get_prosite_raw("PS50948")
record = Prosite.read(handle)
print(record.name)  # Imprime PAN

handle = ExPASy.get_prosite_raw(record.pdoc)
record = Prodoc.read(handle)
