from Bio import Entrez
# Abrir archivo, buscar IDs.

Entrez.email = "dgoretti@lcg.unam.mx"

file = open("file_example.txt", "r")
line = str(file.readlines())
file.close()

line = line.split('"')
saltos = line[2]
saltos = saltos.split(' ')
ids = saltos[4]
ids = ids.replace(']', '')
ids = ids.replace("'", "")
ids = ids.split(",")

if len(ids) >= 3:
    for id in ids:
        # Sacar abstracts
        print('\n' + 'ID: ' + id)
        fetch_handle = Entrez.efetch(db="pubmed", id=id, rettype="abstracts", retmode="text")
        data = fetch_handle.read()  # usar handle.read()
        fetch_handle.close()  # cerrar handle
        print("ABSTRACT: \n")
        print(data)

        results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=id))
        links = results[0]["LinkSetDb"]

        # Sacar citas para cada abstract (for)

        counter = 0
        for link in links:
            print("CITAS: \n")
            if (len(results[0]["LinkSetDb"][counter]["Link"])) > 3:
                print(results[0]["LinkSetDb"][counter]["Link"])
            else:
                print("No hay citas suficientes de este artículo.")
            counter += 1
