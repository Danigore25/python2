'''
NAME
    clases.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com>

DESCRIPTION
    Este programa define una clase principal de insectos y nombra a dos subclases (Lepidoptera y Coleoptera),
    inscribiendo dos organismos, definiendo sus atributos y métodos. Asimismo, se agrega el polimorfismo u overriding
    al método grow.

CATEGORY
    Programación orientada a objetos

USAGE
    clases.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    No se requiere algún paquete.

INPUT
     Datos de los organismos (ya escritos en el código).

OUTPUT
    Datos de los organismos finales en pantalla (de acuerdo con los datos definidos por el código).

EXAMPLES
    Example 1: Se tiene la clase Insect, la cual tiene como atributos generales el ser animal y el ser invertebrado.
    Se define el constructor que será heredado por los organismos pertenecientes a esta superclase, compuesto por las
    características género (taxón), especie (taxón), comida o alimentación (carnívora, herbívora u omnívora) y sexo
    (macho o hembra). Asimismo, se define el método grow mediante el cual se definirá el tamaño del organismo antes de
    crecer y una medida de crecimiento determinada (en centímetros).
    Después se define la primera subclase, Lepidoptera. Esta nueva clase tendrá dos atributos nuevos (el número de alas
    y el ser polinizadores), y reescribirá el cálculo de la longitud derivado del método crecer ya heredado (lo cual
    define a un polimorfismo de tipo overriding). También se define la segunda subclase conocida como Coleoptera, donde
    se va a modificar el método heredado mediante overriding.
    Por último, se postulan dos organismos (una mariposa y una catarina), mencionando sus datos pedidos por herencia a
    partir del constructor de la superclase, así como la longitud del animal pedida por la modificación del método.
    Con los datos del código vemos que imprime:

    {'genre': 'Danaus', 'specie': 'Danaus plexippus', 'food': 'Herbivorous', 'sex': 'Female', 'length': 5.5}
    {'genre': 'Coccinella', 'specie': 'Coccinella septempunctata', 'food': 'Carnivorous', 'sex': 'Male', 'length': 0.6}

    Uso: Este programa podría usarse para saber acerca de organismos como éstos, también podría modificarse el código
    para agregar más especies y más atributos de instancia.

GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/clases.py

'''

# 1. Definir superclase.
class Insect:
    # 1.1. Definir atributos de superclase.
    Animal = True
    Invertebrate = True

    # 1.2. Definir constructor que será heredado.
    def __init__(self, genre, specie, food, sex, length=0):
        self.genre = genre
        self.specie = specie
        self.food = food
        self.sex = sex
        self.length = length

    # 1.3. Definir método.
    def grow(self, last_length):
        self.length += last_length + 5


# 2. Definir subclase Lepidoptera.
class Lepidoptera(Insect):
    # 2.1. Definir atributos de instancia.
    wings = 2
    polinizator = True

    # 2.2. Reescribir método heredado (polimorfismo u overriding).
    def grow(self, last_length):
        self.length += last_length + 0.5


# 3. Definir subclase Coleoptera.
class Coleoptera(Insect):
    # 3.1. Reescribir método heredado (polimorfismo y overriding).
    def grow(self, last_length):
        self.length += last_length + 0.1


# 4. Escribir datos del organismo Butterfly.
Butterfly = Lepidoptera('Danaus', 'Danaus plexippus', 'Herbivorous', 'Female')
Butterfly.grow(5)
# 4.1. Imprimir en pantalla los datos del organismo Butterfly.
print(Butterfly.__dict__)

# 5. Escribir datos del organismo Ladybug.
Ladybug = Coleoptera('Coccinella', 'Coccinella septempunctata', 'Carnivorous', 'Male')
Ladybug.grow(0.5)
# 5.1. Imprimir en pantalla los datos del organismo Ladybug.
print(Ladybug.__dict__)
