import numpy as np
from scipy.linalg import solve

# Données
L = 10  # m
F = 40000  # N
EI = 4e9  # Nm^2

# Système d'équations sous forme Ax = B
A = np.array([[0, 1, 0, 0],  # Équation 1
    [0, 0, L, 1],  # Équation 2
    [L/2, 0, -L/2, -1],  # Équation 3
    [1, 0, -1, 0]    # Équation 4
])

B = np.array([0,
    - (F / (2 * EI)) * (5 * L**2 - L**3 / 6),
    -(F / (12 * EI)) * (L**3 / 8) + (F / (2 * EI)) * (5 * L**2 / 4 - (L / 2)**3 / 6),
    -(F / (12 * EI)) * (3 * (L / 2)**2 ) + (F / (2 * EI)) * (10 * (L / 2) - 3 * (L / 2)**2 / 6)
])

# Résolution du système
solutions = solve(A, B)
C1, C2, C3, C4 = solutions

# Affichage des résultats
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print(f"C3 = {C3}")
print(f"C4 = {C4}")


