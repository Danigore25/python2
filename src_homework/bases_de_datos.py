'''
NAME
    bases_de_datos.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com> <dgoretti@lcg.unam.mx>

DESCRIPTION
    Parte 1. Esta parte del programa pide un correo electrónico para entrar a Entrez y busca los nombres "ECNO" y
    "protein_protein_small_genome" en la base de datos "protein" (dentro de los campos FieldList y LinkList de DbInfo).

    Parte 2. Esta parte del programa pide al usuario el nombre de una persona, palabras que ayuden a buscar uno o más de
    sus investigaciones y una base de datos donde buscar, y se guardan los IDS de los resultados de la búsqueda en el
    archivo ids_author.txt.

CATEGORY
    Búsqueda en bases de datos

USAGE
    bases_de_datos.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    from Bio import Entrez (ambas partes)

INPUT
    Parte 1. Correo electrónico del usuario (por si hay algún problema).

    Parte 2. Nombre del autor(a), dos palabras clave del o los títulos de las obras que se quieren buscar, base de
    datos donde se va a realizar la búsqueda.

OUTPUT
    Parte 1. Datos de la posición, nombre y descripción de "ECNO" y "protein_protein_small_genome":
        Datos de ECNO en FieldList de la base de datos protein:
        Posición:  16 , nombre:  ECNO , descripción:  EC number for enzyme or CAS registry number
        Datos de protein_protein_small genome en LinkList de la base de datos protein:
        Posición:  33 , nombre:  protein_protein_small_genome , descripción:  All proteins from this genome

    Parte 2. Archivo ids_author.txt con lista de IDS de artículos de acuerdo al autor y a las palabras que se buscaron
    en el título (y de acuerdo a la base de datos donde se buscaron).

EXAMPLES
    Parte 1: Se pide el correo electrónico del usuario, y éste lo escribe. Se hace el handle con einfo en la base de
    datos "protein", y se lee el handle con Entrez.read. Se recorren los campos que se encuentran en ["DbInfo"]
    ["FieldList"] de la lectura del handle con un for, y si el campo ["Name"] es igual a "ECNO" se imprime un mensaje
    que contiene la posición en el campo, el nombre y la descripción de "ECNO". Para encontrar
    "protein_protein_small_genome" se realiza una iteración dentro de los campos ["DbInfo"]["LinkList"] de la lectura
    del handle, y al hallar "protein_protein_small_genome" se imprime su posición, su nombre y su descripción.
    De esta forma, se imprime en la pantalla:
        PARTE 1. ECNO y protein_protein_small_genome
        Datos de ECNO en FieldList de la base de datos protein:
        Posición:  16 , nombre:  ECNO , descripción:  EC number for enzyme or CAS registry number
        Datos de protein_protein_small genome en LinkList de la base de datos protein:
        Posición:  33 , nombre:  protein_protein_small_genome , descripción:  All proteins from this genome

    Parte 2: Se tiene el archivo ids_author.txt donde se guardarán los IDs de los trabajos publicados tras realizar la
    búsqueda del autor(a) y las palabras clave; este archivo se abre. Después se pide al usuario el nombre del autor o
    autora que desea buscar, en este ejemplo se agrega el nombre Valeria Mateo-Estrada; luego se pide una palabra del
    título relacionada con una de sus obras y se busca Phylogeographical, y al pedir otra palabra se ingresa
    Epidemiology; al pedir la base de datos donde se va a realizar la búsqueda se escribe pubmed. Se modifican los
    campos para que sea posible realizar la búsqueda y se juntan en una variable llamada termino. Se hace el handle
    usando esearch para buscar la variable termino como el campo term en la base de datos pedida en el input y se hace
    un read del handle. Al final se escriben los IDs de los archivos encontrados con el autor y alguna de las opciones
    de las palabras del título dentro del archivo ids_author.txt.
    Dentro de este archivo queda escrito:
        IDs de la busqueda "Valeria Mateo-Estrada[AUTH] AND (Phylogeographical OR Epidemiology)" en pubmed:
            ['34282943', '32611704']


GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/bases_de_datos.py

'''
# AMBAS PARTES: Importar librería Entrez de Bio.
from Bio import Entrez

# PARTE 1: Búsqueda de "ECNO" y "protein_protein_small_genome en la base de datos Protein"
print("PARTE 1. ECNO y protein_protein_small_genome")

# 1. Pedir correo electrónico.
print("Escriba su correo electrónico: ")
Entrez.email = input()

# 2. Hacer handle con einfo en la base de datos Protein, usar Entrez.read.
handle = Entrez.einfo(db="protein")
record = Entrez.read(handle)

# 3. Recorrer los campos FieldList de DbInfo.
i = -1
for field in record["DbInfo"]["FieldList"]:
    i += 1
    # 3.1. Buscar campo con el nombre "ECNO" e imprimir sus datos.
    if field["Name"] == "ECNO":
        print("Datos de ECNO en FieldList de la base de datos protein: ")
        print("Posición: ", i, ", nombre: ", field["Name"], ", descripción: ", field["Description"])

# 4. Recorrer los campos LinkList de DbInfo.
j = -1
for campo in record["DbInfo"]["LinkList"]:
    j += 1
    # 4.1. Buscar campo con el nombre "protein_protein_small_genome" e imprimir sus datos.
    if campo["Name"] == "protein_protein_small_genome":
        print("Datos de protein_protein_small genome en LinkList de la base de datos protein: ")
        print("Posición: ", j, ", nombre: ", campo["Name"], ", descripción: ", campo["Description"])

# 5. Cerrar handle.
handle.close()


# PARTE 2: Búsqueda de autor(a) y palabras clave del título de sus obras.
print("")
print("PARTE 2. Búsqueda de autor(a) y palabras del título de sus obras")

# 1. Abrir archivo nuevo (se encuentra en el directorio docs).
new = open("../docs/ids_author.txt", "w")

# 2. Pedir nombre del autor o autora.
print("Escriba el nombre de un autor o autora que desee buscar: ")
author = str(input())
author = author + "[AUTH]"

# 3. Pedir la primera palabra del título que se desea buscar.
print("Escriba una palabra del título de una obra de la persona anterior: ")
word1 = str(input())
word1 = " AND (" + word1

# 4. Pedir segunda palabra del título que se desea buscar.
print("Escriba una segunda palabra del título de una obra de la persona anterior: ")
word2 = str(input())
word2 = " OR " + word2 + ")"

# 5. Pedir base de datos.
print("Escriba la base de datos donde desea aplicar su búsqueda: ")
db = str(input())

# 6. Definir campo term de esearch.
termino = author + word1 + word2

# 7. Hacer búsqueda con handle del campo term en la base de datos input. Hacer Entrez.read.
handle = Entrez.esearch(db=db, term=termino)
result = Entrez.read(handle)

# 8. Guardar lista de resultados de búsqueda en archivo.
new.write("IDs de la busqueda \"" + termino + "\" en " + db + ": ")
new.write(str(result["IdList"]))

# 9. Cerrar handle.
handle.close()
