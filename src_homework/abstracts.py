'''
NAME
    abstracts.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com> <dgoretti@lcg.unam.mx>

DESCRIPTION
    Este programa toma los datos de un archivo de texto llamado ids_author.txt (obtenido a través del programa
    bases_de_datos.py), obtiene los IDs que se relacionan con las palabras en el título de un autor y una base de datos
    específica y obtiene el abstract y las citas relacionadas con los IDs.

CATEGORY
    Búsqueda en bases de datos

USAGE
    abstracts.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    from Bio import Entrez

INPUT
    Archivo ids_author.txt previamente llenado a través del archivo bases_de_datos.py

OUTPUT
    ID, abstract de ID y citas de los IDs.

EXAMPLES

    Parte 2: Se tiene el archivo ids_author.txt con la siguiente información (obtenida tras correr el programa
    bases_de_datos.py):

        IDs de la busqueda "Valeria Mateo-Estrada[AUTH] AND (Epidemiology OR bacteria)" en pubmed:
        34282943,32611704,31406982,30625167

    Primero se entra a Entrez por medio del correo electrónico otorgado, se abre el archivo y se leen sus líneas. Se
    obtienen los IDS de la búsqueda por medio de la modificación del texto dentro del archivo. Si hay más de tres IDs se
    imprimen sus abstracts por medio de Entrez.efetch. Se busca con Entrez.read y Entrez.elink las citas, especificando
    el campo LinkName como "pubmed_pmc_refs", con un contador se van recorriendo los links y se evalúa el número de
    citas (si es menor que 3 entonces se imprime que no hay suficientes citas, si es mayor o igual que 3 entonces se
    imprimen las mismas).

    En total hay cuatro IDs, para simplificar el encabezado se describirá lo que se imprime para el ID 31406982:

        ID: 31406982

        ABSTRACT DE ID 31406982:


        1. Genome Biol Evol. 2019 Sep 1;11(9):2531-2541. doi: 10.1093/gbe/evz178.

        Phylogenomics Reveals Clear Cases of Misclassification and Genus-Wide
        Phylogenetic Markers for Acinetobacter.

        Mateo-Estrada V(1), Graña-Miraglia L(1), López-Leal G(1), Castillo-Ramírez S(1).

        Author information:
        (1)Programa de Genómica Evolutiva, Centro de Ciencias Genómicas, Universidad
        Nacional Autónoma de México, Cuernavaca, México.

        The Gram-negative Acinetobacter genus has several species of clear medical
        relevance. Many fully sequenced genomes belonging to the genus have been
        published in recent years; however, there has not been a recent attempt to infer
        the evolutionary history of Acinetobacter with that vast amount of information.
        Here, through a phylogenomic approach, we established the most up-to-date view of
        the evolutionary relationships within this genus and highlighted several cases of
        poor classification, especially for the very closely related species within the
        Acinetobacter calcoaceticus-Acinetobacter baumannii complex (Acb complex).
        Furthermore, we determined appropriate phylogenetic markers for this genus and
        showed that concatenation of the top 13 gives a very decent reflection of the
        evolutionary relationships for the genus Acinetobacter. The intersection between
        our top markers and previously defined universal markers is very small. In
        general, our study shows that, although there seems to be hardly any universal
        markers, bespoke phylogenomic approaches can be used to infer the phylogeny of
        different bacterial genera. We expect that ad hoc phylogenomic approaches will be
        the standard in the years to come and will provide enough information to resolve
        intricate evolutionary relationships like those observed in the Acb complex.

        © The Author(s) 2019. Published by Oxford University Press on behalf of the
        Society for Molecular Biology and Evolution.

        DOI: 10.1093/gbe/evz178
        PMCID: PMC6740150
        PMID: 31406982  [Indexed for MEDLINE]


        CITAS:

        [{'Id': '8437295'}, {'Id': '8409315'}, {'Id': '8407383'}, {'Id': '8269215'}, {'Id': '8226933'},
        {'Id': '8131573'}, {'Id': '8045196'}, {'Id': '7804863'}, {'Id': '7593844'}, {'Id': '7477861'},
        {'Id': '6997731'}]


    NOTA: En dado caso de que el número de IDs sea menor o igual que 3 entonces se imprime lo siguiente:
        No hay suficientes IDs en el archivo ids_author.txt. Intente modificar la búsqueda en el archivo
        "bases_de_datos.py"


GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/abstracts.py

'''

# 1. Importar librería.
from Bio import Entrez

# 2. Poner correo electrónico.

Entrez.email = "dgoretti@lcg.unam.mx"

# 3. Abrir archivo ids_author.txt, obtener líneas.

file = open("../docs/ids_author.txt", "r")
lines = str(file.readlines())
file.close()

# 4. Obtener IDs del archivo de búsqueda.

lines = lines.split('"')
saltos = lines[2]
saltos = saltos.split(' ')
ids = saltos[4]
ids = ids.replace(']', '')
ids = ids.replace("'", "")
ids = ids.split(",")

# 5. Imprimir abstracts de IDs si hay igual o más de tres IDs.

if len(ids) >= 3:
    for id in ids:
        print('\n' + 'ID: ' + id + "\n")

        # 5.1. Usar efetch y leerlo.
        fetch_handle = Entrez.efetch(db="pubmed", id=id, rettype="abstracts", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()
        print("ABSTRACT DE ID " + id + ": \n")

        # 5.2. Imprimir abstracts.
        print(data)

        # 6. Sacar citas.

        # 6.1. Usar read y elink con Entrez para buscar las citas.
        results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=id))
        links = results[0]["LinkSetDb"]

        # 6.2. Imprimir y evaluar si hay citas suficientes.
        counter = 0
        for link in links:
            print("CITAS: \n")
            if (len(results[0]["LinkSetDb"][counter]["Link"])) >= 3:
                print(results[0]["LinkSetDb"][counter]["Link"])
            else:
                print("No hay citas suficientes de este artículo.")
            counter += 1

# 7. Imprimir mensaje si no hay suficientes IDs.
else:
    print("No hay suficientes IDs en el archivo ids_author.txt. Intente modificar la búsqueda en el archivo "
          "bases_de_datos.py")
