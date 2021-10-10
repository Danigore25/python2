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
new = open("file_example.txt", "w")

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
new.write('IDs de la busqueda "' + termino + '" en ' + db + ': ' + '\n')

i = 0
for num in result["IdList"]:
    i += 1
    if i < len(result["IdList"]):
        new.write(str(num) + ',')
    else:
        new.write(str(num))

# 9. Cerrar handle.
handle.close()
