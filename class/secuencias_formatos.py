from Bio.Seq import Seq
from Bio.SeqUtils import nt_search, GC, molecular_weight
from Bio import SeqIO

seqobj = Seq("ATCGATATATACGCGAT")
print(seqobj.translate(to_stop=True))

patron = Seq("ACG")
resultado = nt_search(str(seqobj), patron)
print(resultado)

print(GC(seqobj))
print(molecular_weight(seqobj))

# Ejercicio 1: ORFs
# 1. Se define la secuencia.
sequence = Seq("AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG")

# 2. Se busca el codón de inicio.
inicio = Seq("ATG")
posicion = nt_search(str(sequence), inicio)

# 3. Se recorre la secuencia en busca de los codones de inicio. Se obtiene su secuencia.
for i in range(1, len(posicion)):
    seq_prot = sequence[i:]
    protein = seq_prot.translate(to_stop=True)

# CLASE 2 VERSIÓN 2

# 1. Guardar IDs de records en lista.
mala_calidad = []
umbral = 32
new = open("../docs/records.txt", "w")

# 2. Acceder a Phred_qualities, ir sacando el promedio de cada record.
for record in SeqIO.parse("../docs/sample.fastq", "fastq"):
    promedio = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
    # 2.1. Añadir ID de record si el promedio es menor al umbral de calidad.
    if promedio < umbral:
        mala_calidad.append((promedio, record.id))
    # 2.2. Guardar records que sí superan el umbral.
    if promedio > umbral:
        new.write(record.id)


# 3. Imprimir la longitud de la lista de mala_calidad.
print(len(mala_calidad))

new.close()

# Gen Bank
'''
for gb_record in SeqIO.parse("../docs/aichi.gb", "genbank"):
    print('ID', gb_record.id)
    print('Secuencia', str(gb_record.seq)[0:30], '...')
    print('Longitud', len(gb_record))

    for annotation, value in gb_record.annotations.items():
        print(annotation, value)
'''

# Ejercicio 4.
for gb_record in SeqIO.parse("../docs/virus.gb", "genbank"):
    for annotation, value in gb_record.annotations.items():
        print(annotation, value)

    print(gb_record.annotations['organism'])
    print(gb_record.annotations['sequence_version'])
    print(gb_record.features[0].location)

    # Ejercicio 5. Extraer país y fuente del aislado (source = 0 en features).
    print(gb_record.features[0].qualifiers['isolation_source'])
    print(gb_record.features[0].qualifiers['country'])

    # Guardar inicio y final de la secuencia.
    start = gb_record.features[1].location.nofuzzy_start
    end = gb_record.features[1].location.nofuzzy_end

    # Guardar secuencia dentro del inicio y el final.
    nueva_seq = gb_record.seq[start:end]

    # Traducir proteína.
    protein = nueva_seq.translate()
    print(protein)

    # Imprimir datos del gen L.
    print(gb_record.features[9].qualifiers['gene'])
    start_L = gb_record.features[9].location.nofuzzy_start
    end_L = gb_record.features[9].location.nofuzzy_end
    sequence_L = gb_record.seq[start_L:end_L]
    print(sequence_L)
    rna_L = sequence_L.transcribe()
    print(rna_L[0:5])
    protein_L = sequence_L.translate()
    print(protein_L)

    number = len(gb_record.features) - 1
    while number > -1:
        if gb_record.features[number].qualifiers['gene'] == ['L']:
            print(gb_record.features[number].qualifiers['gene'])
            break
        number -= 1
