import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Constants
Q1 = -6e-6  # Charge on m1 (C)
Q2 = 5e-6  # Charge on m2 (C)
M1 = 0.5   # Mass of m1 (kg)
M2 = 0.6   # Mass of m2 (kg)
L1 = 1.5  # Free length of spring 1 (m)
L2 = 1.6  # Free length of spring 2 (m)
K1 = 50.0  # Spring constant for spring 1 (N/m)
K2 = 50.0  # Spring constant for spring 2 (N/m)
D = 1.0    # Distance between the tops of the springs (m)
G = 9.81  # Acceleration due to gravity (m/s^2)
K_E = 8.9875517923e9  # Coulomb constant (N·m²/C²)

def equilibrium(vars):
    x1, y1, x2, y2 = vars

    r1 = np.sqrt(x1 ** 2 + y1 ** 2)  # Length of spring 1
    r2 = np.sqrt((x2 - D) ** 2 + y2 ** 2)  # Length of spring 2
    r12 = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)  # Distance between charges

    # Spring forces
    F_spring1 = K1 * (r1 - L1)
    F_spring2 = K2 * (r2 - L2)

    # Components of spring forces
    F_spring1_x = F_spring1 * (x1 / r1)
    F_spring1_y = F_spring1 * (y1 / r1)
    F_spring2_x = F_spring2 * ((x2 - D) / r2)
    F_spring2_y = F_spring2 * (y2 / r2)

    # Electrostatic force
    F_electro = K_E * Q1 * Q2 / r12 ** 2
    F_electro_x = F_electro * ((x2 - x1) / r12)
    F_electro_y = F_electro * ((y2 - y1) / r12)

    # Net forces
    eq1 = F_spring1_x + F_electro_x
    eq2 = F_spring1_y - M1 * G + F_electro_y
    eq3 = F_spring2_x - F_electro_x
    eq4 = F_spring2_y - M2 * G - F_electro_y

    return [eq1, eq2, eq3, eq4]

initial_guess = [0.5, 0.5, 1.5, 0.5]

solution = fsolve(equilibrium, initial_guess)

x1, y1, x2, y2 = solution

plt.figure(figsize=(8, 6))
plt.plot([0, x1], [0, -y1], 'b-', label="Spring 1")
plt.plot([D, x2], [0, -y2], 'r-', label="Spring 2")
plt.plot(x1, -y1, 'bo', label="Mass 1")
plt.plot(x2, -y2, 'ro', label="Mass 2")
plt.plot(0, 0, 'ko', label="Anchor 1")
plt.plot(D, 0, 'ko', label="Anchor 2")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Equilibrium Configuration of Masses and Springs")
plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
plt.axvline(0, color='black', linewidth=0.5, linestyle='--')
plt.legend()
plt.grid()
plt.text(x1, -y1, f"({x1:.2f}, {y1:.2f})", fontsize=10, color='blue', ha='left')
plt.text(x2, -y2, f"({x2:.2f}, {y2:.2f})", fontsize=10, color='red', ha='left')
plt.show()
