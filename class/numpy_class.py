import numpy as np

# Array de primera dimensión
array_1D = np.array([1, 2, 3])
print(array_1D)

# Array bidimensional
array_2D = np.array([[1, 2, 3], (2, 3, 4)])
print(array_2D)

# Array tridimensional
array_3D = np.array([[[1, 2], [3, 4]],   [[5, 6], [7, 8]]])
print(array_3D)

# Array de una dimensión de unidades de absorbancia
ecoli_matraz = np.array([0.1, 0.15, 0.19, 0.5,
                         0.9, 1.4, 1.8, 2.1, 2.3])
print(ecoli_matraz.ndim)
print(ecoli_matraz.shape)
print(len(ecoli_matraz))

# Biomasa en unidades de absorbancia  (OD600)
ecoli_m_b = np.array([[0.1, 0.15, 0.19, 0.5,  # Matraz 250 mL
                       0.9, 1.4, 1.8, 2.1, 2.3],
                      [0.1, 0.17, 0.2, 0.53,  # Biorreactor 50 L
                       0.97, 1.43, 1.8, 2.1,  2.8],
                      [0.1, 0.17, 0.2, 0.52,  # B. alimentado 50 L
                       0.95, 1.41, 1.8, 2.2,  2.8]])
print(ecoli_m_b.ndim)
print(ecoli_m_b.shape)
print(len(ecoli_m_b))

# Si una bacteria equivale a 0.39:
print(ecoli_m_b*0.39)

# Operaciones
produccion = np.array([[16, 14], [12, 9]])
print(produccion+produccion)  # Suma de 2 biorreactores (por bacteria y metabolito)
print(np.sum(produccion*2, axis=0))  # Suma total por metabolito
print(np.sum(produccion*2, axis=1))  # Por bacteria

# Si se contaminó la mitad del líquido de uno de los reactores:
total = produccion*2
contaminado = produccion/2
total_real = total - contaminado
print(total_real)

# EJERCICIO 1 ---------------------------------------------
# produccion = np.array([[5, 3], [11, 7], [4, 9], [2, 6]])
# Sacar la transpuesta para ir multiplicando cada producción
# t_produccion = produccion.T
# print(t_produccion)
# costos = np.array([3.5, 5, 7, 4.3])
# print(costos/t_produccion)

# TIPOS DE DATOS ----------------------------------------------------------------------------------------------
from sys import getsizeof

np_float = np.array([1.0, 2.0, 3.0, 4.0])
print("Tipo de dato\t", np_float.dtype, "\nTamaño en bytes\t", getsizeof(np_float))

np_int = np.array([1, 2, 3, 4])
print("Tipo de dato\t", np_int.dtype, "\nTamaño en bytes\t", getsizeof(np_int))

# Array booleano
bool_np = np.array([True, False, True, False])
print(bool_np.dtype)
# Acceder con array booleano
print(np_int[bool_np])

# costo_unitario = (costos/t_produccion).T
# menor_costo = costo_unitario.min()
# mayor_costo = costo_unitario.max()

# menor_costo_bol = costo_unitario == menor_costo
# mayor_costo_bol = costo_unitario == mayor_costo
# Acceder a ambos:
# print(costo_unitario[(mayor_costo_bol | menor_costo_bol)])

# Arrays estructurados
mascotas = np.array([('Freya', 6, 6.5), ('Senna', 1, 2.5)], dtype=[('nombre', (np.str_, 10)), ('edad', np.int32),
                                                                   ('peso', np.float64)])
print(mascotas)

sort_age = np.sort(mascotas, order='edad')
sort_name = np.sort(mascotas, order='nombre')

print(sort_age)
print(sort_name)


# 1. Importar librería numpy
import numpy as np


# 2. Generar array de producción
produccion = np.array([(5, 3), (11, 7), (4, 9), (2, 6)], dtype=[('30 ° C', np.int32), ('35 ° C', np.int32)])

print(produccion)
# 3. Transponer matriz de producción
t_produccion = produccion.T

print(t_produccion)
# 4. Generar array de costos
costos = np.array([3.5, 5, 7, 4.3], dtype=[('Costo', np.int32)])

print(costos)
# 5. Generar array de costo unitario
costo_unitario = costos/t_produccion

print(costo_unitario)
# costo_unitario_array = np.array(costo_unitario, dtype=[('30 °C', np.str_), ('35 ° C', np.str_)])
