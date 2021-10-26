import pandas as pd
import numpy as np

serie = pd.Series(['a', 'b', 'c', 'd', 'e'],
                  index=['a', 'b', 'c', 'd', 'e'],
                  name="Ejemplo Serie")

print(serie)


ecoli_matraz = pd.Series([0.1, 0.15, 0.19, 0.5,
                         0.9, 1.4, 1.8, 2.1, 2.3],
                         index=['t1', 't2', 't3', 't4',
                                't5', 't6', 't7', 't8', 't9'],
                         name='Matraz')
print(ecoli_matraz)

ODs = pd.Series([0.2, 0.2, 0.4, 0.1, 0.2, 0.1, 0.2, 0.4, 0.1],
                index=[8, 4, 1, 2, 3, 0, 5, 7, 6],
                name='Ajustes')

# EJERCICIO 1 ----------------------------------------------------------------------

produccion = pd.Series([5, 11, 4, 7, 2], index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'])

costos = pd.Series([5, 4.3, 7, 3.5], index=['gen1', 'gen2', 'gen3', 'gen5'])

costo_unitario = costos/produccion.T

print(costo_unitario)
print(costo_unitario.min())

# -----------------------------------------------------
nan_test = pd.Series([0.1, None, 2.1, 2.3], name='Matraz')
print(nan_test.count())

# loc y iloc
series_test = pd.Series([5.1, 2.2, 1.1, 3.1, 4.2], index=[5, 2, 1, 3, 4])
print(series_test)
print(series_test.loc[1])
print(series_test.iloc[1])

# EJERCICIO 2 ------------------------------------------------------------------
bool_min = costo_unitario == costo_unitario.min()
bool_max = costo_unitario == costo_unitario.max()

print(costo_unitario[bool_min | bool_max])

# Repetir índices
regulon = pd.Series(['aidB', 'alaS', 'accB', 'accC', 'bhsA'], index=['AidB', 'AlaS', 'AccB', 'AccB', 'ComR'],
                    name='Genes regulados')
print(regulon.loc['AccB'])
print(regulon.loc['AidB'])


# Clases en series
class Mamifero:
    vertebrado = True

    def haz_ruido(self):
        print('aaaaaaaaaaaaaaaaaaaaaaaaaaa')


array_clase = pd.Series([np.sum, 'a', Mamifero], name='objetos')
jerbo = array_clase.iloc[2]
print(jerbo.haz_ruido())
