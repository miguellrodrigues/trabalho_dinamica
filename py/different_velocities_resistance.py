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

# dados foguete
m = .875  # massa do foguete
r = .15  # raio do bico (m)

# Arrasto
A = np.pi * r**2
cd = 2.327  # coeficiente de arrasto
rho = 1.1839
g = 9.80665  # aceleracao gravitacional

V_0 = 120  # modulo da velocidade inicial
velocities = np.linspace(1, 100, 1000)

theta = np.pi/4

V0_X = velocities * np.cos(theta)  # componente horizontal da velocidade inicial
V0_Y = velocities * np.sin(theta)  # componente vertical da velociade inicial

theorical_horizontal_ranges = 2*V0_X*m/(A*cd*rho) - 2*V0_X*m*np.exp(-(A*V0_Y*cd*rho + 2*g*m*lambertw(-np.sqrt((A**2*V0_Y**2*cd**2*rho**2 + 4*A*V0_Y*cd*g*m*rho + 4*g**2*m**2)*np.exp(-A*V0_Y*cd*rho/(g*m)))*np.exp(-1)/(2*g*m)) + 2*g*m)/(2*g*m))/(A*cd*rho)
theorical_vertical_ranges = 2*m*(A*V0_Y*cd*rho - 2*g*m*np.log(A*V0_Y*cd*rho + 2*g*m) + np.log(2**(2*g*m)*g**(2*g*m)*m**(2*g*m)))/(A**2*cd**2*rho**2)

real_part_thr = np.real(theorical_horizontal_ranges)

ax, fig1 = plt.subplots()
fig1.set_title('Alcances horizontais')
fig1.set_xlabel('Velocidade inicial (m/s)')
fig1.set_ylabel('Alcance (m)')
fig1.plot(velocities, real_part_thr, label='Alcance horizontal')
plt.savefig('images/air_resistance/alcances_horizontais_velocidade.png', dpi=300)

ax, fig2 = plt.subplots()
fig2.set_title('Alcances verticais')
fig2.set_xlabel('Velocidade inicial (m/s)')
fig2.set_ylabel('Alcance (m)')
fig2.plot(velocities, theorical_vertical_ranges, label='Alcance vertical')

plt.savefig('images/air_resistance/alcances_verticais_velocidade.png', dpi=300)
plt.show()
