lista = []
print("Introduzca el número de genes que desea buscar: ")
numero = int(input())

while numero > 0:
    print("Introduzca el gen que desee guardar en la lista: ")
    element = str(input())
    lista.append(element)
    numero -= 1

print("Los valores en la lista son: ", lista)
