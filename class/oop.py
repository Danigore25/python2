class Animal:
    name = 'Ninguno'
    years = 'None'

    def __init__(self, alimentacion, peso, altura):
        self.alimentacion = alimentacion
        self.peso = peso
        self.altura = altura

    def hazruido(self, ruido):
        print(ruido)


class Perro(Animal):
    def __init__(self, alimentacion, peso, altura, hazruido):
        Animal.__init__(self, alimentacion, peso, altura)
        self.hazruido = hazruido


class Gato(Animal):
    usa_arenero = True

    def __init__(self, alimentacion, peso, altura, hazruido):
        Animal.__init__(self, alimentacion, peso, altura)
        self.hazruido = hazruido


chihuahua = Perro('omnivoro', '2', '30', 'guau')
chihuahua.hazruido()

angora = Gato('carnivoro', '5', '100', 'miau')
angora.hazruido()
