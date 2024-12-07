import numpy as np
from scipy import linalg
import random


A1 = np.array([
    [2.554219275, 0.871733993, 0.052575899, 0.240740262, 0.316022841],
    [0.871733993, 0.553460938, -0.070921727, 0.255463951, 0.707334556],
    [0.052575899,-0.070921727, 3.409888776, 0.293510439, 0.847758171],
    [0.240740262, 0.255463951, 0.293510439, 1.108336850, -0.206925123],
    [0.316022841, 0.707334556, 0.847758171, -0.206925123, 2.374094162],
])
A2 = np.array([
    [2.645152285, 0.544589368, 0.009976745, 0.327869824, 0.424193304],
    [0.544589368, 1.730410927, 0.082334875, -0.057997220, 0.318175706],
    [0.009976745, 0.082334875, 3.429845092, 0.252693077, 0.797083832],
    [0.327869824, -0.057997220, 0.252693077, 1.191822050, -0.103279098],
    [0.424193304, 0.318175706, 0.797083832, -0.103279098, 2.502769647]
])


b = [-0.642912346, -1.408195475, 4.595622394, -5.073473196, 2.178020609]


result1 = linalg.solve(A1, b)
result2 = linalg.solve(A2, b)

war_wl_A1, wec_wl_A1 = np.linalg.eig(A1)
war_wl_A2, wec_wl_A2 = np.linalg.eig(A2)

print("Wynik A1y=b dla b: ")
print(result1)
print("\nWynik A2y=b dla b: ")
print(result2)
print("\nWartosci wlasne A1: " + str(war_wl_A1))
print("\nWectory wlasne A1: " + str(wec_wl_A1))
print("\nWartosci wlasne A2: " + str(war_wl_A2))
print("\nWektory wlasne A2: "+ str(wec_wl_A2))

rows = 100
cols = 5

matrix_delta_b = np.zeros((rows, cols))

for i in range(rows):
    for j in range(cols):
        matrix_delta_b[i, j] = random.uniform(9.6e-7, 10.4e-7)


for i in range(rows):
    d = np.zeros((1, cols))
    delta_b = matrix_delta_b[i, :]
    d = np.add(b, delta_b)
    result_zab1 = linalg.solve(A1, d)
    result_zab2 = linalg.solve(A2, d)

    print("\nWynik A1y = b + delta(b): ")
    print(result_zab1)
    print("\nWynik A2y = b + delta(b): ")
    print(result_zab2)
