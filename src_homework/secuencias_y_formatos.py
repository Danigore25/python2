'''
NAME
    secuencias_y_formatos.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com>

DESCRIPTION
    Este programa recibe un archivo GenBank como input, le pide al usuario los genes que quiere buscar dentro del mismo
    e imprime el nombre del organismo, la versión de la secuencia, la fuente del aislado, el país, y para cada gen
    muestra su nombre, los primeros 15 nucleótidos de DNA y RNA que tiene y sus primeros 15 aminoácidos.

CATEGORY
    Archivos de GenBank

USAGE
    secuencias_y_formatos.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    Se requiere importar la librería SeqIO de Bio.

INPUT
     Archivo GenBank (en este caso se utiliza el archivo virus.gb). Asimismo, se pide el número y los nombres de los
     genes que se desean buscar.

OUTPUT
    Datos generales del organismo:
        Organismo:  <Nombre del organismo>
        Versión de la secuencia:  <Versión>
        Fuente del aislado:  <Fuente>
        País:  <País de origen>

    Para cada uno de los genes:
        Nombre del gen:  <Gen que se busca>
        Primeros quince nucleótidos de DNA:  <Secuencia de DNA>
        Primeros quince nucleótidos de RNA:  <Secuencia de RNA>
        Primeros quince aminoácidos:  <Secuencia proteica>

EXAMPLES
    Example 1: Se tiene el archivo ../docs/virus.gb. Se escribe la ruta de este archivo de GenBank al ser pedido como
    input. En este caso se quieren buscar 2 genes, por lo que se escribe ese valor cuando se pida el número de genes
    que se desean buscar; en este caso, se nombran dos de los genes dentro del archivo virus.gb: N y L, estos nombres
    se agregan a una lista que los contendrá. Se tiene una función llamada resumen que tomará la lista de los genes que
    se van a buscar y el archivo en donde queremos buscarlos. Aquí se nombra al registro como gb_record llamando a
    SeqIO.parse para meternos dentro de los datos del archivo de GenBank. Después se van imprimiendo cada uno de los
    datos mencionados (el nombre del organismo, la versión de la secuencia, la fuente del aislado y el país, y para
    imprimir los datos de los genes éstos se buscan dentro del archivo y se imprime su nombre, sus secuencias de DNA,
    RNA y proteínas resumidas a 15 caracteres). Para ir buscando los genes se tiene que realizar otro loop, donde se
    vaya buscando cada gen en la lista features de gb_record. Para esto se utiliza una variable llamada maximo, la cual
    va a ayudar a buscar estos genes por medio del índice indicado en la lista features. Al final, se llama a la
    función resumen.
    En este caso, los resultados son:
        Organismo:  Isfahan virus
        Versión de la secuencia:  1
        Fuente del aislado:  ['Phlebotomus papatasi']
        País:  ['Iran:Isfahan province']
        Nombre del gen:  ['N']
        Primeros quince nucleótidos de DNA:  ATGACTTCTGTAGTA
        Primeros quince nucleótidos de RNA:  AUGACUUCUGUAGUA
        Primeros quince aminoácidos:  MTSVVKRIATGSSVL
        Nombre del gen:  ['L']
        Primeros quince nucleótidos de DNA:  ATGGATGAGTACTCT
        Primeros quince nucleótidos de RNA:  AUGGAUGAGUACUCU
        Primeros quince aminoácidos:  MDEYSEEKWGDSDEE


GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/secuencias_y_formatos.py

'''

# 1. Importar librerías.
from Bio import SeqIO

# 2. Pedir ruta del archivo GenBank.
print("Escriba la ruta del archivo: ")
archivo_input = input()

# 3. Pedir número de genes que se van a buscar.
genes_input = []
print("Introduzca el número de genes que desea buscar: ")
number = int(input())

# 4. Pedir los nombres de los genes e incluirlos en una lista.
while number > 0:
    print("Introduzca el nombre gen que desea buscar: ")
    element = str(input())
    genes_input.append(element)
    number -= 1


# 5. Realizar función resumen.
def resumen(archivo, genes):
    """
    :param archivo: Archivo de GenBank pedido por el usuario.
    :param genes: Lista de genes que se buscarán en el archivo.
    :return: Datos del archivo de GenBank y de los genes que se buscaron en él.
    """
    # 5.1. Hacer loop que pida los datos del archivo llamando a SeqIO.parse.
    for gb_record in SeqIO.parse(archivo, "genbank"):
        print("Organismo: ", gb_record.annotations['organism'])
        print("Versión de la secuencia: ", gb_record.annotations['sequence_version'])
        print("Fuente del aislado: ", gb_record.features[0].qualifiers['isolation_source'])
        print("País: ", gb_record.features[0].qualifiers['country'])

        # 5.2. Hacer loop que busque los genes en la lista dentro del archivo.
        for gene in genes:
            # 5. 2. 1. Definir una variable que recorra la lista gb_record.features.
            maximo = len(gb_record.features) - 1
            # 5. 2. 2. Hacer loop de búsqueda de los genes para imprimir sus datos.
            while maximo > -1:
                if gb_record.features[maximo].qualifiers['gene'] == [gene]:
                    print("Nombre del gen: ", gb_record.features[maximo].qualifiers['gene'])
                    start = gb_record.features[maximo].location.nofuzzy_start
                    end = gb_record.features[maximo].location.nofuzzy_end
                    dna = gb_record.seq[start:end]
                    print("Primeros quince nucleótidos de DNA: ", dna[0:15])
                    rna = dna.transcribe()
                    print("Primeros quince nucleótidos de RNA: ", rna[0:15])
                    protein = dna.translate()
                    print("Primeros quince aminoácidos: ", protein[0:15])
                    break
                maximo -= 1


# 6. Llamar a la función resumen con los datos pedidos por el usuario.
resumen(archivo_input, genes_input)
