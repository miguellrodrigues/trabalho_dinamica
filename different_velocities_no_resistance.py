import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import lambertw

plt.style.use(
    [
        'science',
        'nature',
        'grid',
    ]
)

b = .125
m = 1  # massa do foguete
g = 9.80665  # aceleracao gravitacional

Velocities = np.linspace(10, 100, 1000)
theta = np.pi/4

X_0 = .0  # posicao horizontal inicial
Y_0 = 1e-3 # posicao vertical inicial

V0_X = Velocities * np.cos(theta)
V0_Y = Velocities * np.sin(theta)

theorical_horizontal_ranges = (V0_X / g) * (V0_Y + np.sqrt(V0_Y**2 + 2*g*Y_0))
theorical_vertical_ranges = Y_0 + (V0_Y**2 / (2*g))

ax, fig1 = plt.subplots()
fig1.set_title('Alcances horizontais')
fig1.set_xlabel('Velocidade inicial (m/s)')
fig1.set_ylabel('Alcance (m)')
fig1.plot(Velocities, theorical_horizontal_ranges, label='Alcance horizontal')
plt.savefig('images/no_air_resistance/alcances_horizontais_velocidade.png', dpi=300)

ax, fig2 = plt.subplots()
fig2.set_title('Alcances verticais')
fig2.set_xlabel('Velocidade inicial (m/s)')
fig2.set_ylabel('Alcance (m)')
fig2.plot(Velocities, theorical_vertical_ranges, label='Alcance vertical')

plt.savefig('images/no_air_resistance/alcances_verticais_velocidade.png', dpi=300)
plt.show()
