'''
NAME
    numpy.py

VERSION
    1.0

AUTHOR
    Daniela Goretti Castillo León <danigore22@gmail.com>

DESCRIPTION
    Este programa crea arrays estructurados a partir de arrays vistos en la clase de Biopython.

CATEGORY
    Creación de arrays

USAGE
    numpy.py [sin opciones]

ARGUMENTS
    No se requieren argumentos.

PACKAGES
    import numpy as np

INPUT
    Arrays dados por el código (costos, produccion).

OUTPUT
    Impresiones de los arrays en forma estructurada.

EXAMPLES
    Example 1: Se tiene el array de producción [[5, 3], [11, 7], [4, 9], [2, 6]], así como el array de costos por gen
    [3.5, 5, 7, 4.3]; se escriben ambos arrays de manera estructurada. Para sacar el array de costo unitario se debe
    transponer el array de producción, para después dividir el array de costos entre la misma. Al final, se obtiene un
    array de dos dimensiones, el cual va a transformarse en un array estructurado.

GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/numpy.py

'''

# 1. Importar librería numpy
import numpy as np


# 2. Generar array de producción
produccion = np.array([(5, 3), (11, 7), (4, 9), (2, 6)], dtype=[('30 ° C', np.int32), ('35 ° C', np.int32)])
print(produccion)

# 3. Transponer matriz de producción
t_produccion = produccion.T

# 4. Generar array de costos
costos = np.array([3.5, 5, 7, 4.3], dtype=[('Costo', np.int32)])
print(costos)

# 5. Generar array de costo unitario
costo_unitario = (costos/t_produccion).T
costo_unitario_array = np.array(costo_unitario, dtype=[('30 °C', np.int32), ('35 ° C', np.int32)])
print(costo_unitario_array)
