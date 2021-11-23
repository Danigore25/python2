import pandas as pd
import numpy as np

pd_DF = pd.DataFrame(np.random.rand(3, 2),
                     columns=["columna_1", "columna_2"],
                     index=['a', 'b', 'c'])

# A partir de objeto Series:
produccion = pd.Series([5, 11, 4, 7, 2],
                       index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'],
                       name='production')
costos = pd.Series([5, 4.3, 7, 3.5], index=['gen1', 'gen2', 'gen3', 'gen5'], name='costos')

costo_beneficio = pd.DataFrame({'costos': costos, 'produccion': produccion})

# Con loc y iloc:
print(costo_beneficio.loc['gen1'::2, 'costos'])

costo_beneficio['doble'] = costo_beneficio.costos + costo_beneficio.produccion
print(costo_beneficio)

# EJERCICIO 1 -----------------------------------------------------------------------------------------
produccion = pd.Series([5, 11, 4, 7, 2], index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'], name='produccion')
costos = pd.Series([3.5, 5, 7, 4.3], index=['gen1', 'gen2', 'gen3', 'gen5'], name='costos')
costo_beneficio2 = pd.DataFrame({'costos': costos, 'produccion': produccion})
print(costo_beneficio2)

# Columna costo unitario
costo_beneficio2['unitario'] = costo_beneficio2.costos/costo_beneficio2.produccion
print(costo_beneficio2)

# Obtener los índices de valores máximos de cada columna
print(costo_beneficio2.idxmax())

# EJERCICIO 2 -----------------------------------------------------------------------------------------------

produccion_30 = pd.Series([5, 11, 4, 7, 2], index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'], name='produccion')
produccion_35 = pd.Series([3, 7, 9, 4, 6], index=['gen1', 'gen2', 'gen3', 'gen4', 'gen5'], name='produccion')
costos = pd.Series([3.5, 5, 7, 4.3], index=['gen1', 'gen2', 'gen3', 'gen5'], name='costos')
costo_beneficio3 = pd.DataFrame({'costos': costos, 'produccion 30°': produccion_30, 'produccion 35°': produccion_35})
print("Costo beneficio 3 \n", costo_beneficio3)

# Obtener costo unitario en una sola operación
columnas_interes = ['produccion 30°', 'produccion 35°']
producciones = costo_beneficio3.loc[:, columnas_interes]

# Division subset entre costos
costos_unitarios = producciones.div(costo_beneficio3.costos, axis=0)

# Renombrar
costos_unitarios.rename(columns={'produccion 30°': 'costo unitario 30°', 'produccion 35°': 'costo unitario 35°'},
                        inplace=True)

# Añadir columna organismos
organismos = np.random.choice(['procariotas', 'eucariotas', 'arqueas'], 5, p=[0.5, 0.3, 0.2])
costo_beneficio3['organismos'] = organismos
print(costo_beneficio3)

