# Librerias a usar
from Bio import ExPASy, Entrez, SwissProt
from Bio.ExPASy import Prosite, Prodoc

# Correo electronico

Entrez.email = "dgoretti@lcg.unam.mx"

# Crear una función que tome una lista de terminos GO y una lista de IDs de UniProt


def get_swissprot(go, ids):
    # Abrir archivo donde se guardara la informacion

    # Extraer archivos SwissProt de cada ID
    for ID in ids:
        handle = ExPASy.get_sprot_raw(ID)
        record = SwissProt.read(handle)
        print(record.__dict__.keys())

        # Buscar si algun termino GO aparece en el archivo
        for GO in go:
            # En caso de encontrarlo
            if record.cross_references[2][1] == GO:
                # Obtener la informacion
                print('El ID es: ' + ID + '\n')
                print('Nombre de proteina: ' + record.description.split(';')[0] + '\n')
                print('Nombre y definicion del GO: ' + GO + record.cross_references[2][2] + '\n')
                print('Organismo: ' + record.organism + '\n')

                # FALTA: ----------------------------------------------------------------
                # Imprimir localización subcelular:
                for comment in record.comments:
                    type_c = comment.split(':')
                    if type_c[0] == 'SUBCELLULAR LOCATION':
                        print('Localización subcelular: ', comment, '\n')

                # Abstract de fuente - Sacar ID de Pubmed
                # print('Abstract de una de las fuentes: ')
                # fetch_handle = Entrez.efetch(db="pubmed", id= id_pubmed, rettype="abstract", retmode="text")
                # data = fetch_handle.read()
                # fetch_handle.close()
                # print(data)

                # Imprimir PROSITE
                for reference in record.cross_references:
                    if 'PROSITE' in reference:
                        print(reference[1])

                        # handle = ExPASy.get_prosite_raw(reference[1])
                        # record = Prosite.read(handle)
                        # print(record.name)  # Imprime PAN

                        # handle = ExPASy.get_prosite_raw(record.pdoc)
                        # record = Prodoc.read(handle)
                        # print(record.text)


GO_Terms = ["GO:0044178", "GO:0046755", "GO:0046761", "GO:0046760", "GO:0039702", "GO:0046765", "GO:0046762"]
uniprot_IDs = ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC", "POLG_YEFV1", "POLG_DEN1W", "Q6W352_CVEN9", "D9SV67_CLOC7",
               "A9KSF7_LACP7", "B8I7R6_RUMCH"]

# Pedir y crear archivo nuevo
# print("Escriba el nombre o ruta del archivo donde desea escribir (sin terminación): ")
# file_input = str(input() + ".txt")

get_swissprot(GO_Terms, uniprot_IDs)
