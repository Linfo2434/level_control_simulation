import numpy as np
import matplotlib.pyplot as plt

# Intervalo de tiempo y particion
a = 0
b = 100
T = (b * 1000)

h = (b - a) / T
t = np.linspace(a, b, T + 1)

x_prop = np.zeros((2, T + 1))
v_prop = np.zeros(T + 1)
Fin_prop = np.zeros(T + 1)
Fout_prop = np.zeros(T + 1)

x_PI = np.zeros((2, T + 1))
v_PI = np.zeros(T + 1)
Fin_PI = np.zeros(T + 1)
Fout_PI = np.zeros(T + 1)

x_ST = np.zeros((2, T + 1))
v_ST = np.zeros(T + 1)
Fin_ST = np.zeros(T + 1)
Fout_ST = np.zeros(T + 1)

x_TP = np.zeros((2, T + 1))
v_TP = np.zeros(T + 1)
Fin_TP = np.zeros(T + 1)
Fout_TP = np.zeros(T + 1)

x_MD = np.zeros((2, T + 1))
v_MD = np.zeros(T + 1)
Fin_MD = np.zeros(T + 1)
Fout_MD = np.zeros(T + 1)


# Parámetros del sistema
AT = 1
beta = 1000
Cv = 10
g0 = 9.8
alpha = 0.08
k = 0.1
k_MD = 0.008
kp = 0.1
ki = 0.1
k_TP = 3 # k para tiempo predefinido
xsp = 8.1


# Condiciones iniciales
x_prop[0,0],x_prop[1,0] = 1,0
Fin_prop[0] = 2.27 

x_PI[0,0],x_PI[1,0] = 1,0
Fin_PI[0] = 2.27 

x_ST[0,0],x_ST[1,0] = 1,0
Fin_ST[0] = 2.27 

x_TP[0,0],x_TP[1,0] = 1,0
Fin_TP[0] = 2.27 

x_MD[0,0],x_MD[1,0] = 1,0
Fin_MD[0] = 2.27 


# Trayectoria de referencia
yd = np.array([0, 0])

# Acción de control Proporcional
u = lambda t, x: k * (x - xsp)
def f1(t, x, v, Fin):
    return np.array([(Fin - Cv * v * (1.45 * 10 ** (-4) * beta * g0 * x[0]) ** (0.5)) / AT,
                     alpha * np.sign(x[0] - xsp)])

# Runge Kutta de 4to orden (Proporcional)
for i in range(len(t) - 1):
    Fin_prop[i + 1] = np.random.uniform(2.27, 6.8) if t[i] % 15 == 0 else Fin_prop[i]
    v_prop[i] = 0 if (u(t[i], x_prop[0, i]) + x_prop[1, i]) < 0 or (u(t[i], x_prop[0, i]) + x_prop[1, i]) > 1 else (u(t[i], x_prop[0, i]) + x_prop[1, i])

    k_1 = f1(t[i], x_prop[:, i], v_prop[i], Fin_prop[i])
    k_2 = f1(t[i] + 0.5 * h, x_prop[:, i] + 0.5 * h * k_1, v_prop[i], Fin_prop[i])
    k_3 = f1(t[i] + 0.5 * h, x_prop[:, i] + 0.5 * h * k_2, v_prop[i], Fin_prop[i])
    k_4 = f1(t[i] + h, x_prop[:, i] + k_3 * h, v_prop[i], Fin_prop[i])
    x_prop[:, i + 1] = x_prop[:, i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h

    Fout_prop[i] = Cv * v_prop[i] * (1.45 * 10 ** (-4) * beta * g0 * x_prop[0, i]) ** (0.5) / AT

# Acción de control PI
u_PI = lambda t, x: kp * abs(x - xsp) + ki * (x - xsp) * t
def f2(t, x, v, Fin):
    return np.array([(Fin - Cv * v * (1.45 * 10 ** (-4) * beta * g0 * x[0]) ** (0.5)) / AT,
                     alpha * np.sign(x[0] - xsp)])
# Runge Kutta de 4to orden (PI)
for i in range(len(t) - 1):
    Fin_PI[i + 1] = np.random.uniform(2.27, 6.8) if t[i] % 15 == 0 else Fin_PI[i]
    v_PI[i] = 0 if (u_PI(t[i], x_PI[0, i]) + x_PI[1, i]) < 0 or (u_PI(t[i], x_PI[0, i]) + x_PI[1, i]) > 1 else (u_PI(t[i], x_PI[0, i]) + x_PI[1, i])

    k_1 = f2(t[i], x_PI[:, i], v_PI[i], Fin_PI[i])
    k_2 = f2(t[i] + 0.5 * h, x_PI[:, i] + 0.5 * h * k_1, v_PI[i], Fin_PI[i])
    k_3 = f2(t[i] + 0.5 * h, x_PI[:, i] + 0.5 * h * k_2, v_PI[i], Fin_PI[i])
    k_4 = f2(t[i] + h, x_PI[:, i] + k_3 * h, v_PI[i], Fin_PI[i])
    x_PI[:, i + 1] = x_PI[:, i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h

    Fout_PI[i] = Cv * v_PI[i] * (1.45 * 10 ** (-4) * beta * g0 * x_PI[0, i]) ** (0.5) / AT

# Accion de control Super Twisting
u_ST = lambda t, x: k * abs(x - xsp) ** (0.5) * np.sign(x - xsp)
def f3(t, x, v, Fin):
    return np.array([(Fin - Cv * v * (1.45 * 10 ** (-4) * beta * g0 * x[0]) ** (0.5)) / AT,
                     alpha * np.sign(x[0] - xsp)])

# Runge Kutta de 4to orden (ST)
for i in range(len(t) - 1):
    Fin_ST[i + 1] = np.random.uniform(2.27, 6.8) if t[i] % 15 == 0 else Fin_ST[i]
    v_ST[i] = 0 if (u_ST(t[i], x_ST[0, i]) + x_ST[1, i]) < 0 or (u_ST(t[i], x_ST[0, i]) + x_ST[1, i]) > 1 else (u_ST(t[i], x_ST[0, i]) + x_ST[1, i])

    k_1 = f3(t[i], x_ST[:, i], v_ST[i], Fin_ST[i])
    k_2 = f3(t[i] + 0.5 * h, x_ST[:, i] + 0.5 * h * k_1, v_ST[i], Fin_ST[i])
    k_3 = f3(t[i] + 0.5 * h, x_ST[:, i] + 0.5 * h * k_2, v_ST[i], Fin_ST[i])
    k_4 = f3(t[i] + h, x_ST[:, i] + k_3 * h, v_ST[i], Fin_ST[i])
    x_ST[:, i + 1] = x_ST[:, i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h

    Fout_ST[i] = Cv * v_ST[i] * (1.45 * 10 ** (-4) * beta * g0 * x_ST[0, i]) ** (0.5) / AT

# Accion de control Tiempo Predefinido
u_TP = lambda t, x: -(k_TP * (x - xsp)) / (t - b)
def f4(t, x, v, Fin):
    return np.array([(Fin - Cv * v * (1.45 * 10 ** (-4) * beta * g0 * x[0]) ** (0.5)) / AT,
                     alpha * np.sign(x[0] - xsp)])

# Runge Kutta de 4to orden (TP)
for i in range(len(t) - 1):
    Fin_TP[i + 1] = np.random.uniform(2.27, 6.8) if t[i] % 15 == 0 else Fin_TP[i]
    v_TP[i] = 0 if (u_TP(t[i], x_TP[0, i]) + x_TP[1, i]) < 0 or (u_TP(t[i], x_TP[0, i]) + x_TP[1, i]) > 1 else (u_TP(t[i], x_TP[0, i]) + x_TP[1, i])

    k_1 = f4(t[i], x_TP[:, i], v_TP[i], Fin_TP[i])
    k_2 = f4(t[i] + 0.5 * h, x_TP[:, i] + 0.5 * h * k_1, v_TP[i], Fin_TP[i])
    k_3 = f4(t[i] + 0.5 * h, x_TP[:, i] + 0.5 * h * k_2, v_TP[i], Fin_TP[i])
    k_4 = f4(t[i] + h, x_TP[:, i] + k_3 * h, v_TP[i], Fin_TP[i])
    x_TP[:, i + 1] = x_TP[:, i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h

    Fout_TP[i] = Cv * v_TP[i] * (1.45 * 10 ** (-4) * beta * g0 * x_TP[0, i]) ** (0.5) / AT

# Accion de control Modos deslizantes
u_MD = lambda t, x: k_MD * np.sign(x - xsp)
def f5(t, x, v, Fin):
    return np.array([(Fin - Cv * v * (1.45 * 10 ** (-4) * beta * g0 * x[0]) ** (0.5)) / AT,
                     alpha * np.sign(x[0] - xsp)])
# Runge Kutta de 4to orden (MD)
for i in range(len(t) - 1):
    Fin_MD[i + 1] = np.random.uniform(2.27, 6.8) if t[i] % 15 == 0 else Fin_MD[i]
    v_MD[i] = 0 if (u_MD(t[i], x_MD[0, i]) + x_MD[1, i]) < 0 or (u_MD(t[i], x_MD[0, i]) + x_MD[1, i]) > 1 else (u_MD(t[i], x_MD[0, i]) + x_MD[1, i])

    k_1 = f5(t[i], x_MD[:, i], v_MD[i], Fin_MD[i])
    k_2 = f5(t[i] + 0.5 * h, x_MD[:, i] + 0.5 * h * k_1, v_MD[i], Fin_MD[i])
    k_3 = f5(t[i] + 0.5 * h, x_MD[:, i] + 0.5 * h * k_2, v_MD[i], Fin_MD[i])
    k_4 = f5(t[i] + h, x_MD[:, i] + k_3 * h, v_MD[i], Fin_MD[i])
    x_MD[:, i + 1] = x_MD[:, i] + (1 / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h

    Fout_MD[i] = Cv * v_MD[i] * (1.45 * 10 ** (-4) * beta * g0 * x_MD[0, i]) ** (0.5) / AT

# Graficas
#Grafica Proporcional
fig = plt.figure('Proporcional',figsize=(12,8))


ax1 = plt.subplot(2,2,1)
ax1.plot(t,x_prop[0,:],label = "$h(t)$")
ax1.grid()
ax1.legend()
ax1.set_ylim([8.09,8.11])

ax2 = plt.subplot(2,2,3)
ax2.plot(t,x_prop[1,:],label = "$x_2(t)$")
ax2.grid()
ax2.legend()

ax3 = plt.subplot(2,2,2)
ax3.plot(t[:-1],v_prop[:-1],label="u_TP(t)")
ax3.grid()
ax3.legend()

ax4 = plt.subplot(2,2,4)
ax4.plot(t,Fin_prop,label="$F_{in}(t)$")
ax4.plot(t[:-1],Fout_prop[:-1],label="$F_{out}(t)$",alpha = 0.7)
ax4.set_ylim([0,10])
ax4.grid()
ax4.legend()

plt.show()


#Gráfica PI
fig = plt.figure('P.I.',figsize=(12,8))


ax1 = plt.subplot(2,2,1)
ax1.plot(t,x_PI[0,:],label = "$h(t)$")
ax1.grid()
ax1.legend()
ax1.set_ylim([8.09,8.11])

ax2 = plt.subplot(2,2,3)
ax2.plot(t,x_PI[1,:],label = "$x_2(t)$")
ax2.grid()
ax2.legend()

ax3 = plt.subplot(2,2,2)
ax3.plot(t[:-1],v_PI[:-1],label="u(t)")
ax3.grid()
ax3.legend()

ax4 = plt.subplot(2,2,4)
ax4.plot(t,Fin_PI,label="$F_{in}(t)$")
ax4.plot(t[:-1],Fout_PI[:-1],label="$F_{out}(t)$",alpha = 0.7)
ax4.set_ylim([0,10])
ax4.grid()
ax4.legend()

plt.show()



#Grafica TP
fig = plt.figure('Control de Tiempo Predefinido',figsize=(12,8))


ax1 = plt.subplot(2,2,1)
ax1.plot(t,x_TP[0,:],label = "$h(t)$")
ax1.grid()
ax1.legend()
ax1.set_ylim([8.09,8.11])

ax2 = plt.subplot(2,2,3)
ax2.plot(t,x_TP[1,:],label = "$x_2(t)$")
ax2.grid()
ax2.legend()

ax3 = plt.subplot(2,2,2)
ax3.plot(t[:-1],v_TP[:-1],label="u(t)")
ax3.grid()
ax3.legend()

ax4 = plt.subplot(2,2,4)
ax4.plot(t,Fin_TP,label="$F_{in}(t)$")
ax4.plot(t[:-1],Fout_TP[:-1],label="$F_{out}(t)$",alpha = 0.7)
ax4.set_ylim([0,10])
ax4.grid()
ax4.legend()

plt.show()

#Grafica ST
fig = plt.figure('Super Twisting',figsize=(12,8))


ax1 = plt.subplot(2,2,1)
ax1.plot(t,x_ST[0,:],label = "$h(t)$")
ax1.grid()
ax1.legend()
ax1.set_ylim([8.09,8.11])

ax2 = plt.subplot(2,2,3)
ax2.plot(t,x_ST[1,:],label = "$x_2(t)$")
ax2.grid()
ax2.legend()

ax3 = plt.subplot(2,2,2)
ax3.plot(t[:-1],v_ST[:-1],label="u_ST(t)")
ax3.grid()
ax3.legend()

ax4 = plt.subplot(2,2,4)
ax4.plot(t,Fin_ST,label="$F_{in}(t)$")
ax4.plot(t[:-1],Fout_ST[:-1],label="$F_{out}(t)$",alpha = 0.7)
ax4.set_ylim([0,10])
ax4.grid()
ax4.legend()

plt.show()

#Gráfica MD
fig = plt.figure('Modos deslizantes',figsize=(12,8))


ax1 = plt.subplot(2,2,1)
ax1.plot(t,x_MD[0,:],label = "$h(t)$")
ax1.grid()
ax1.legend()
ax1.set_ylim([8.09,8.11])

ax2 = plt.subplot(2,2,3)
ax2.plot(t,x_MD[1,:],label = "$x_2(t)$")
ax2.grid()
ax2.legend()

ax3 = plt.subplot(2,2,2)
ax3.plot(t[:-1],v_MD[:-1],label="u(t)")
ax3.grid()
ax3.legend()

ax4 = plt.subplot(2,2,4)
ax4.plot(t,Fin_MD,label="$F_{in}(t)$")
ax4.plot(t[:-1],Fout_MD[:-1],label="$F_{out}(t)$",alpha = 0.7)
ax4.set_ylim([0,10])
ax4.grid()
ax4.legend()

plt.show()
