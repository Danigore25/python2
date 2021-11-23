from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons


def dySIS(y, t, lamda, mu):  # SI/SIS model
    dy_dt = lamda*y*(1-y)-mu*y
    return dy_dt


def dySIR(y, t, lamda, mu):  # SIR model,
    i, s = y
    di_dt = lamda*s*i-mu*i
    ds_dt = -lamda*s*i
    return np.array([di_dt, ds_dt])

# Parametrización
number = 1e5  # total number of people
lamda = 0.2  # Daily contact rate, the average number of susceptible persons who are effectively in contact with the sick each day
sigma = 2.5  # Number of contacts during infectious period
mu = lamda/sigma  # Daily cure rate, the ratio of the number of patients cured each day to the total number of patients
tEnd = 200  # Forecast date length
t = np.arange(0.0, tEnd, 1)  # (start,stop,step)
i0 = 1e-4  # Initial value of the proportion of patients
s0 = 1-i0  # Initial value of the proportion of susceptible persons
Y0 = (i0, s0)  # Initial value of the differential equation system

# Resolver ecuaciones diferenciales
ySI = odeint(dySIS, i0, t, args=(lamda, 0))  # SI model
ySIS = odeint(dySIS, i0, t, args=(lamda, mu))  # SIS model
ySIR = odeint(dySIR, Y0, t, args=(lamda, mu))  # SIR model

# Gráfica
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.4)
plt.title("Comparison among SI, SIS and SIR models")
plt.xlabel('time')
plt.axis([0, tEnd, -0.1, 1.1])
si_plt, = plt.plot(t, ySI, ':g', label='i(t)-SI')
sis_plt, = plt.plot(t, ySIS, '--g', label='i(t)-SIS')
sir_i_plt, = plt.plot(t, ySIR[:, 0], '-r', label='i(t)-SIR')
sir_s_plt, = plt.plot(t, ySIR[:, 1], '-b', label='s(t)-SIR')
sir_r_plt, = plt.plot(t, 1-ySIR[:, 0]-ySIR[:, 1], '-m', label='r(t)-SIR')
plt.legend(loc='best')

plt.show()

# GRAFICA INTERACTIVA ---
# Agregamos las barras interactivas
axcolor = 'lightgoldenrodyellow'

# Generamos el área de las barras
axlambda = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axsigma = plt.axes([0.25, 0.18, 0.65, 0.03], facecolor=axcolor)
axi0 = plt.axes([0.25, 0.26, 0.65, 0.03], facecolor=axcolor)

# Agregamos la información
slambda = Slider(axlambda, 'Daily contact rate', 0.1, 1,
               valinit=lamda, color="green")
ssigma = Slider(axsigma, 'Contacts during\ninfectious period', 0.1, 10,
                valinit=sigma)
si0 = Slider(axi0, 'Initial proportion\nof patients', 1e-4, 5e-1,
                valinit=i0, color="orange")

# Conectar información

def update(val, ):
    lamda = slambda.val
    sigma = ssigma.val
    i0 = si0.val
    mu = lamda / sigma
    s0 = 1 - i0
    Y0 = (i0, s0)
    ySI = odeint(dySIS, i0, t, args=(lamda, 0))  # SI model
    ySIS = odeint(dySIS, i0, t, args=(lamda, mu))  # SIS model
    ySIR = odeint(dySIR, Y0, t, args=(lamda, mu))  # SIR model
    si_plt.set_ydata(ySI)
    sis_plt.set_ydata(ySIS)
    sir_i_plt.set_ydata(ySIR[:, 0])
    sir_s_plt.set_ydata(ySIR[:, 1])
    sir_r_plt.set_ydata(1 - ySIR[:, 0] - ySIR[:, 1])
    fig.canvas.draw_idle()
plt.show()

slambda.on_changed(update)
ssigma.on_changed(update)
si0.on_changed(update)

# Botones para ver solo un modelo de
rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('SI', 'SIS', 'SIR'), active=0)
lines = {'SI': [si_plt], 'SIS': [sis_plt],
         'SIR': [sir_i_plt, sir_s_plt, sir_r_plt]}


def select_model(label):
    # La línea seleccionada no es transparente
    for line_m in lines[label]:
        line_m.set_alpha(1)

    # Las demás líneas serán transparentes
    for others in set(lines.keys()) - set([label]):
        for line_m in lines[others]:
            line_m.set_alpha(0)
    fig.canvas.draw_idle()


radio.on_clicked(select_model)

plt.show()

