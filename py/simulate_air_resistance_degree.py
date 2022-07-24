from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


plt.style.use(
    [
        'nature',
        'science',
        'notebook',
    ]
)


r = .15  # raio do bico (m)
h = .1  # altura do bico (m)

# Arrasto
A = np.pi * r**2
cd = .327  # coeficiente de arrasto
rho = 1.1839

b = .5 * rho * cd * A

m = .875  # massa do foguete
g = 9.80665  # aceleracao gravitacional

V_0 = 200  # modulo da velocidade inicial (m/s)
thetas = np.arange(0, 90, 1)  # angulo de lancamento (rad)

simulation_time = 60  # duracao da simulacao
simulation_step = 1e-3  # tempo de simulacao

iterations = int(simulation_time / simulation_step)

time = np.linspace(.0, simulation_time, iterations)
time_span = [time[0], time[-1]]


def dSdt(_, S):
    _, vx, y, vy = S

    if y < 0:
        vy = 0

    return [
        vx,
        ((-b * vx**2) / m),
        vy,
        ((-b * vy**2) / m) - g
    ]


horizontal_ranges = np.zeros((len(thetas)))
vertical_ranges   = np.zeros((len(thetas)))

for i in range(len(thetas)):
    theta = np.deg2rad(thetas[i])

    V0_X = V_0 * np.cos(theta)
    V0_Y = V_0 * np.sin(theta)

    sol = solve_ivp(
        dSdt,
        time_span,
        [1e-3, V0_X, 1e-3, V0_Y],
        t_eval=time
    )

    x = sol.y[0]
    y = sol.y[2]

    y_zero_index = np.argmin(np.abs(y))
    time_to_reach_zero = time[y_zero_index]

    # alcance horizontal
    horizontal_range = x[y_zero_index]

    # alcance vertical
    max_y_index = np.argmax(y[:y_zero_index])
    vertical_range = y[max_y_index]

    horizontal_ranges[i] = horizontal_range
    vertical_ranges[i] = vertical_range


# polyfit for horizontal range
p_h = np.polyfit(thetas, horizontal_ranges, 32)
p_v = np.polyfit(thetas, vertical_ranges, 6)

# evaluate the polynomial
y_h = np.polyval(p_h, thetas)
y_v = np.polyval(p_v, thetas)

np.set_printoptions(precision=6, suppress=True)

print("\nHorizontal range polynomial:")
print(p_h)
print("\nVertical range polynomial:")
print(p_v)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

ax.set_xlabel("Ângulo de lançamento (graus)")
ax.grid(True)

ax.plot(thetas, horizontal_ranges, 'o', label='Alcance horizontai', color='red')
ax.set_ylabel('Alcance horizontal (m)')

ax1 = ax.twinx()
ax1.plot(thetas, vertical_ranges, 'o', label='Alcance vertical', color='green')
ax1.grid(False)
ax1.set_ylabel('Alcance vertical (m)')
plt.legend(loc='lower center')
plt.savefig('sim_data.png', dpi=300)

fig1, ax_1 = plt.subplots(1, 1, figsize=(8, 6))

ax_1.grid(True)

ax_1.set_xlabel("Ângulo de lançamento (graus)")
ax_1.plot(thetas, y_h, '-', color='red', label='Aproximação polinomial (Horizontal)')
ax_1.set_ylabel('Alcance horizontal (m)')

ax11 = ax_1.twinx()
ax11.grid(False)

ax11.plot(thetas, y_v, '-', label='Aproximação polinomial (Vertical)', color='green')
ax11.set_ylabel('Alcance vertical (m)')
plt.legend(loc='lower center')

plt.savefig('poly_approx.png', dpi=300)
plt.show()
