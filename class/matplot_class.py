import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


# creamos figura
fig = plt.figure()
# creamos un eje
ax = plt.axes()

x = np.linspace(0, 10, 1000)
# ploteamos (x, y respectivamente)
ax.plot(x, np.sin(x))

ax.set(xlim=(0, 10), ylim=(-2, 2),  # limites
       xlabel="x", ylabel="sen(x)",  # etiquetas
       title="grafiquita")       # titulo

# al momento de usar un plt show , matplotlib descarta la figura,
plt.show()

# subplots
# aqui creamos dos, distribuidos(1,2)

fig, axs = plt.subplots(nrows=1, ncols=2)

axs[0].plot(x, np.sin(x))  # primer axes

axs[1].plot(x, np.cos(x))  # segundo axes

# poner labels, etc

# enseña plot
plt.show()

# plt.plot con marcadores
plt.plot(x, np.sin(x), "-.")
plt.plot(x, np.cos(x), "o")
plt.show()
