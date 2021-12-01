import argparse
import Bio.Entrez
# import numpy as np
# import matplotlib.pyplot as plt

Bio.Entrez.email = "dgoretti@lcg.unam.mx"

print("Introduce los terminos de busqueda")
print("Ingrese el pais en el que desea buscar")
print("Ingrese el año de publicacion que le interesa")

# Paso de parametros por argumentos
parser = argparse.ArgumentParser(description="Busqueda de terminos en pubmed")
parser.add_argument("-w", "--word",
                    metavar="word of interest",
                    type=str,
                    help="Terminos a buscar",
                    required=True)

# Search by country and year of publication.
parser.add_argument("-c", "--country",
                    metavar="country of publication",
                    type=str,
                    help="Pais donde se desean buscar las publicaciones",
                    required=True)

parser.add_argument("-y", "--year",
                    metavar="year of publication",
                    type=int,
                    help="Fecha de publicacion del articulo",
                    required=True)

args = parser.parse_args()

search2 = args.country + "[CNTY] AND " + args.year + "[PDAT]"

# through handle the terms in "search" will be looked for in pubmed database
handle = Bio.Entrez.esearch(db="pubmed", term=args.word)

# The results that coincide with the terms in "search" are stored in result
result = Bio.Entrez.read(handle)

# in archiv, the list with the corresponding IDs will be stored
archiv = result["IdList"]
handle.close()
