import argparse
from Bio import Entrez
Entrez.email = "dgoretti@lcg.unam.mx"

print("Introduce los terminos de busqueda")

# Paso de parametros por argumentos
parser = argparse.ArgumentParser(description="Busqueda de terminos en pubmed")
parser.add_argument("-i", "--input",
                    metavar="word of interest",
                    type=str,
                    help="Terminos a buscar",
                    required=True)

args = parser.parse_args()
search = args.input


# through handle the terms in "search" will be looked for in pubmed database
handle = Entrez.esearch(db="pubmed", term=search)

# The results that coincide with the terms in "search" are stored in result
result = Entrez.read(handle)

# in archiv, the list with the corresponding IDs will be stored
archiv = result["IdList"]
handle.close()
