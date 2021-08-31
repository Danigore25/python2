class animal():
    name = 'Ninguno'
    years = 'None'

    def __init__(self, alimentacion, peso, altura):
        self.alimentacion = alimentacion
        self.peso = peso
        self.altura = altura

    def hazruido(self, fuerte=True):
        if fuerte:
            print('AAAAAAAAH')
        else:
            print('aaaaah')


class perro(animal):
    ruido = 'ladra'


class gato(animal):
    accion = 'usa_arenero'
