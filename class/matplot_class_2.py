from sklearn.datasets import load_iris
import matplotlib.pyplot as plt

# cargar set de datos
iris = load_iris()
# transponemos datos de iris
caracteristicas = iris.data.T

# PLOOOOT!!!!
plt.scatter(caracteristicas[0], caracteristicas[1], alpha=0.2,
            s=100*caracteristicas[3], c=iris.target, cmap="viridis")

# etiquetas de plot
plt.title("scatter plot exploratorio")
plt.xlabel(iris.feature_names[0])
plt.ylabel(iris.feature_names[1])

plt.show()

# SLIDER ----------------------------------------------------------------------
import numpy as np
from matplotlib.widgets import Slider
from math import pi

# parametros de slider
a = 5  # eje x
b = 5  # eje y
# crear datos
t = np.linspace(0, 2*pi, 100)
# elipse centro 0,0
x = a * np.cos(t)
y = b * np.sin(t)

# creando figure y axes
fig, ax = plt.subplots()
# ajustar los espacio de figure para slider
plt.subplots_adjust(bottom=0.3)
# plotear datos
p, = ax.plot(x, y, color="red")

# espacio para slider (x,y ancho, altura)
recuadro = plt.axes([0.25, 0.1, 0.65, 0.03])
recuadro2 = plt.axes([0.25, 0.2, 0.65, 0.03])

# crear slider
factor = Slider(recuadro, 'factor a ', valmin=0.1, valmax=10, valinit=1, valstep=1)
factor2 = Slider(recuadro2, 'factor b', valmin=0.1, valmax=10, valinit=8, valstep=1)

# funcion para actualizar datos cada que cambie slider


def update(val):
    # 1. tomamos valor de slider
    nuevo_valor = factor.val
    nuevo_valor2 = factor2.val

    # 2. actualizamos operaciones que son afectadas por el factor
    x = nuevo_valor * np.cos(t)
    y = nuevo_valor2 * np.sin(t)

    # 3. Modificamos datos que están en plot (p.)
    p.set_xdata(x)
    p.set_ydata(y)

    # 4. dibujar linea con actualizaciones
    plt.draw()

# llamar a update cada que cambie el valor de slider 1


factor.on_changed(update)

# llamar a update cada que cambie el valor de slider 2
factor2.on_changed(update)

plt.show()

# RADIO BUTTONS --------------------------------------------------------------------------------
from matplotlib.widgets import RadioButtons


# creando funciones a graficar
t = np.linspace(0,10,1000)
f_seno = np.sin(t)
f_coseno = np.cos(t)
f_tang = np.tan(t)

# creamos figure y axes
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.3)
# plot de datos
plot, = ax.plot(t, f_seno, color="green")
# recuadro
axcolor = 'lightgoldenrodyellow'
rax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)

# crear radio buttons
radio = RadioButtons(rax, ("seno", "coseno", "tangente"))

# funcion que actualiza la grafica


def tipo_funcion(label):
    fun_dict = {"seno": f_seno, "coseno": f_coseno, "tangente": f_tang}
    nueva_funcion = fun_dict[label]
    # cambiar datos de y
    plot.set_ydata(nueva_funcion)
    # plotear nueva funcion
    plt.draw()


# cada que haya cambio en radio button se llama a funcion de actualizacion
radio.on_clicked(tipo_funcion)

plt.show()
