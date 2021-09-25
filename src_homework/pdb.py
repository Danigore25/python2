'''
NAME
    pdb.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com>

DESCRIPTION
    Este programa pregunta al usuario la ruta de un archivo PDB y la cadena en la que desea buscar un aminoácido
    especificado por el mismo. Al final imprime la lista de los residuos que se encontraron (y si la lista está vacía
    imprime más bien un mensaje).

CATEGORY
    Análisis de archivos PDB

USAGE
    pdb.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    from Bio import PDB

INPUT
     Ruta de archivo PDB, cadena en la que se desea buscar y residuo de aminoácido deseado.

OUTPUT
    Lista de aminoácidos encontrados. Si no se encuentra ese aminoácido en la cadena se imprime un mensaje:
    Si se encontraron:
        Se encontraron  <número>  residuos  <aminoácido>  en la cadena  <nombre de la cadena> .
        Lista de los residuos  <aminoácido>  encontrados en la cadena  <nombre de la cadena>  en el archivo input:
        <Lista>

    Si no se encontraron:
        No se encontraron los residuos  <aminoácido>  en la cadena especificada. Intente con otra cadena u otro residuo.

EXAMPLES
    Example 1: Se escribe la ruta del archivo ../docs/1fat.pdb. Vamos a buscar en la cadena A el aminoácido CYS, por lo
    que anotamos ambas cosas cuando se especifique en la pantalla. Estos datos se llaman a partir de la función
    lista_residuos, donde se crea la lista donde se guardarán los residuos de aminoácido que se encuentran en la cadena.
    Se quitan los warnings al archivo PDB y se crea el objeto struct. Se va iterando en los modelos del objeto, después
    en las cadenas y si la cadena es la misma que la del input entonces se sigue buscando. Se va iterando en la cadena
    hasta buscar el residuo de aminoácido que se escogió, y si se encuentra entonces se añade a la lista. En este caso
    este aminoácido no se encuentra, por lo que se imprime:

        No se encontraron los residuos  CYS  en la cadena especificada. Intente con otra cadena u otro residuo.

    Example 2: En la misma ruta del archivo se va a buscar el aminoácido TYR en la cadena A. Estos datos entran a la
    función lista_residuos, se buscan en el archivo PDB y se imprime lo siguiente en pantalla:

        Se encontraron  5  residuos  TYR  en la cadena  A .
        Lista de los residuos  TYR  encontrados en la cadena  A  en el archivo input:
        [<Residue TYR het=  resseq=5 icode= >, <Residue TYR het=  resseq=51 icode= >, <Residue TYR het=  resseq=127
        icode= >, <Residue TYR het=  resseq=167 icode= >, <Residue TYR het=  resseq=180 icode= >]



GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/pdb.py

'''

# 1. Importar librerías
from Bio import PDB

# 2. Pedir ruta del archivo.
print("Escriba la ruta del archivo PDB (ejemplo: ../docs/1fat.pdb): ")
archivo_input = input()

# 3. Pedir nombre de la cadena en el archivo.
print("Escriba la cadena que desea buscar en el archivo (ejemplo: A): ")
cadena_input = input()

# 4. Pedir nombre del residuo de aminoácido.
print("Escriba el nombre del residuo de aminoácido que se va a buscar en la cadena (ejemplo: CYS)")
aminoacido_input = input()


# 5. Crear función para obtener los residuos.
def lista_residuos(archivo, cadena, aminoacido):
    """
    :param archivo: Archivo input PDB.
    :param cadena: Cadena input en la que se van a buscar los residuos.
    :param aminoacido: Residuo de aminoácido input.
    :return: Lista de residuos de aminoácido en la cadena (si se encuentran, si no entonces se imprime un mensaje).
    """
    # 6. Crear lista donde se guardarán los residuos.
    residues = list()

    # 7. Quitar warnings en PDB, hacer objeto struct.
    parser = PDB.PDBParser(QUIET=True)
    obj_struc = parser.get_structure("prot_pdb", archivo)

    # 8. Hacer iteración para guardar los residuos de aminoácidos.
    for model in obj_struc:
        for chain in model:
            if chain.id == cadena:
                for residuo in chain:
                    if residuo.get_resname() == aminoacido:
                        residues.append(residuo)

    # 9. Imprimir lista si la lista de residuos tiene elementos, si no, se imprime otro mensaje.
    if len(residues) == 0:
        print("No se encontraron los residuos ", aminoacido, " en la cadena especificada. Intente con otra cadena u "
                                                             "otro residuo.")

    else:
        print("Se encontraron ", len(residues), " residuos ", aminoacido, " en la cadena ", cadena, ".")
        print("Lista de los residuos ", aminoacido, " encontrados en la cadena ", cadena, " en el archivo input: ")
        print(residues)


# 10. Llamar a función.
lista_residuos(archivo_input, cadena_input, aminoacido_input)
