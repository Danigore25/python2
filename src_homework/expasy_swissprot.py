# Librerias a usar
from Bio import ExPASy
from Bio import SwissProt

# Crear una función que tome una lista de terminos GO y una lista de IDs de UniProt


def get_swissprot(go, ids, my_file):
    # Abrir archivo donde se guardara la informacion
    file = open(my_file, "w")

    # Extraer archivos SwissProt de cada ID
    for ID in ids:
        handle = ExPASy.get_sprot_raw(ID)
        record = SwissProt.read(handle)

        # Buscar si algun termino GO aparece en el archivo
        for GO in go:
            # En caso de encontrarlo
            if record.cross_references[2][1] == GO:
                # Obtener la informacion
                file.write('El ID es: ' + ID + '\n')
                file.write('Nombre de proteina: ' + record.description.split(';')[0] + '\n')
                file.write('Nombre y definicion del GO: ' + GO + record.cross_references[2][2] + '\n')
                file.write('Organismo: ' + record.organism + '\n')


GO_Terms = ["GO:0046755", "GO:0046761", "GO:0046760", "GO:0039702", "GO:0046765", "GO:0046762"]
uniprot_IDs = ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC", "POLG_YEFV1", "POLG_DEN1W", "Q6W352_CVEN9", "D9SV67_CLOC7",
               "A9KSF7_LACP7", "B8I7R6_RUMCH"]

# Pedir y crear archivo nuevo
print("Escriba el nombre o ruta del archivo donde desea escribir (sin terminación): ")
file_input = str(input() + ".txt")

get_swissprot(GO_Terms, uniprot_IDs, file_input)
