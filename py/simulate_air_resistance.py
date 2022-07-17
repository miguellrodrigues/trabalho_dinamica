from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


plt.style.use(
    [
        'nature',
        'science',
        'grid',
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

V_0 = 120  # modulo da velocidade inicial (m/s)
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
p_h = np.polyfit(np.deg2rad(thetas), horizontal_ranges, 32)
p_v = np.polyfit(np.deg2rad(thetas), vertical_ranges, 8)

# evaluate the polynomial
y_h = np.polyval(p_h, np.deg2rad(thetas))
y_v = np.polyval(p_v, np.deg2rad(thetas))

np.set_printoptions(precision=6, suppress=True)

print("\nHorizontal range polynomial:")
print(p_h)
print("\nVertical range polynomial:")
print(p_v)


plt.plot(thetas, horizontal_ranges, 'o', label='horizontal range')
plt.plot(thetas, y_h, '-', label='horizontal range polynomial')
plt.plot(thetas, vertical_ranges, 'o', label='vertical range')
plt.plot(thetas, y_v, '-', label='vertical range polynomial')
plt.legend()
plt.show()
