import numpy as np
import matplotlib.pyplot as plt


def is_approximately(x, y, epsilon):
    return abs(x - y) <= epsilon



# body type: sphere
b = .047  # drag coefficient


m = 1 # mass
g = 9.80665  # gravity acceleration


# initial conditions
v = (45, np.radians(90))  # velocity

vx_0 = v[0] * np.cos(v[1])
vy_0 = v[0] * np.sin(v[1])

x0, y0 = .0, .0  # initial position


def x_t(t):
    return ((-m * vx_0) / b) * np.exp(-(b*t)/m) + x0 + (m/b) * vx_0


def v_xt(t):
    return vx_0 * np.exp(-(b*t)/m)


def y_t(t):
    return ((-m*g)/b) * t - ((m*vy_0)/b) * np.exp(-(b*t)/m) - ((m**2 * g)/(b**2)) * np.exp(-(b*t)/m) + y0 + ((m*vy_0)/b) + ((m**2 * g)/(b**2))


def v_yt(t):
    return ((-m*g)/b) + (vy_0 + (m*g)/b) * np.exp(-(b*t)/m)


simulation_time = 15
simulation_step = 1e-3

iterations = int(simulation_time / simulation_step)
time_values = np.linspace(0, simulation_time, iterations)


x_values = np.zeros(iterations)
y_values = np.zeros(iterations)

x_velocities = np.zeros(iterations)
y_velocities = np.zeros(iterations)

x_values[0] = x0
y_values[0] = y0

x_velocities[0] = vx_0
y_velocities[0] = vy_0

time_stop = 0

for i in range(1, iterations):
    yt = y_t(time_values[i])
    
    if yt <= 1e-2:
        time_stop = i
        break

    x_values[i] = x_t(time_values[i])
    y_values[i] = yt

    x_velocities[i] = v_xt(time_values[i])
    y_velocities[i] = v_yt(time_values[i])


print(' ')
print('Time in air:', time_values[time_stop])

# find the time the ball take to have a velocity of 0
time_to_zero = 0
for i in range(1, time_stop):
    if is_approximately(y_velocities[i], 0, 1e-2):
        time_to_zero = time_values[i]
        break

print(' ')
print('Time to reach max y value:', time_to_zero)
print('Time to reach the ground (from max y value):', time_values[time_stop] - time_to_zero)
print(' ')


fig, axs = plt.subplots(2, 2)

axs[0][0].plot(time_values[:time_stop], x_values[:time_stop], color='blue')
axs[1][0].plot(time_values[:time_stop], y_values[:time_stop], color='red')

axs[0][1].plot(time_values[:time_stop], x_velocities[:time_stop], color='blue')
axs[1][1].plot(time_values[:time_stop], y_velocities[:time_stop], color='red')

axs[0][0].legend(['x(t)'])
axs[1][0].legend(['y(t)'])

axs[0][1].legend(['vx(t)'])
axs[1][1].legend(['vy(t)'])

fig1, axs1 = plt.subplots(1, 1)

axs1.plot(x_values[:time_stop], y_values[:time_stop], color='black')
axs1.legend(['y(x)'])

plt.show()