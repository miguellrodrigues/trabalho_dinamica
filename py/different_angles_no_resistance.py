import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import lambertw

plt.style.use(
    [
        'nature',
        'science',
        'grid',
        'notebook',
    ]
)

m = .875  # massa do foguete
g = 9.80665  # aceleracao gravitacional

V_0 = 120  # modulo da velocidade inicial
thetas = np.linspace(0, 90, 1000)

X_0 = .0  # posicao horizontal inicial
Y_0 = 1e-3 # posicao vertical inicial

V0_X = V_0 * np.cos(np.deg2rad(thetas))  # componente horizontal da velocidade inicial
V0_Y = V_0 * np.sin(np.deg2rad(thetas))  # componente vertical da velociade inicial

theorical_horizontal_ranges = (V0_X / g) * (V0_Y + np.sqrt(V0_Y**2 + 2*g*Y_0))
theorical_vertical_ranges = Y_0 + (V0_Y**2 / (2*g))

ax, fig1 = plt.subplots()
fig1.set_title('Alcances horizontais')
fig1.set_xlabel('Angulo de lançamento (graus)')
fig1.set_ylabel('Alcance (m)')
fig1.plot(thetas, theorical_horizontal_ranges, label='Alcance horizontal')
plt.savefig('images/no_air_resistance/alcances_horizontais_angulo.png', dpi=300)

ax, fig2 = plt.subplots()
fig2.set_title('Alcances verticais')
fig2.set_xlabel('Angulo de lançamento (graus)')
fig2.set_ylabel('Alcance (m)')
fig2.plot(thetas, theorical_vertical_ranges, label='Alcance vertical')

plt.savefig('images/no_air_resistance/alcances_verticais_angulo.png', dpi=300)
plt.show()
