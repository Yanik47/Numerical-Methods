
import time
import numpy as np
from scipy import linalg
from scipy.linalg import det
import mpmath
import matplotlib.pyplot as plt

mpmath.mp.dps = 50
N = 124
epsilon = 10e-8


macierz_bib = np.zeros((N, N))
x_bib = [0.0] * N
for i in range(N):
    for j in range(N):
        if i == j:
            macierz_bib[i][j] = 3.0
        elif i == j + 1 or j == i + 1:
            macierz_bib[i][j] = 1.0
        elif i == j + 2 or j == i + 2:
            macierz_bib[i][j] = 0.15

    x_bib[i] = i+1.0

print("Macierz incjalizowana jako 2-wymiarowa tablica: ")
print(macierz_bib)
print("Wektor x: ")
print(x_bib)
start = time.time()
wyznacznik_strt = det(macierz_bib)
result = linalg.solve(macierz_bib, x_bib)
stop = time.time()
wynik = float(stop - start)
print("Czas, za ktory rozwiazany byl uklad rownan: ")
print(wynik)
print("Wyznacznik liczony przez det(A): ")
print(wyznacznik_strt)
print("Wynik(wektor y) liczony z linalg.solve: ")
print(result)

x_start = [0.0]*N

x_result = []
A = np.zeros((5, N))
D = [0.0]*N
    
x = [0.0]*N
b = [0.0]*N

A[0][0] = 3
A[1][0] = 1
A[2][0] = 0.15

A[0][1] = 1
A[1][1] = 3
A[2][1] = 1
A[3][1] = 0.15

for i in range(2, N-2):
    #Wypelniamy macierz A
    A[0][i] = 0.15
    A[1][i] = 1.0
    A[2][i] = 3.0
    A[3][i] = 1.0
    A[4][i] = 0.15
    #Wektor b
    b[i] = i+1
    #maicewrz diagonalna -D^-1
    D = -1/A[2][i]

A[0][N-2] = 0.15
A[1][N-2] = 1
A[2][N-2] = 3
A[3][N-2] = 1

A[0][N-1] = 0.15
A[1][N-1] = 1
A[2][N-1] = 3



N_macierz = np.zeros((2, N))
M_macierz = np.zeros((3, N))
for i in range(N):
    #liczenie macierzy N oraz M dla Gaussa
    N_macierz[0][i] = -A[0][i]
    N_macierz[1][i] = -A[1][i]
    M_macierz[0][i] = A[2][i]
    M_macierz[1][i] = A[3][i]
    M_macierz[2][i] = A[4][i]
    #Liczenie macierzy M = -D^-1(L + U) dla Jacobiego
for i in range(N):
    A[0][i] = A[0][i]*D[i]
    A[1][i] *= D[i]
    A[2][i] = 0
    A[3][i] *= D[i]
    A[4][i] *= D[i]


def metodaJacobiego(N, x_start, epsilon, result, A, D, b):
    x = [0.0]*N
    x[0] = D[0]*(-1)*b[0] + A[3][0]*x_start[1] + A[4][0]*x_start[2]
    x[1] = D[1]*(-1)*b[1] + A[0][1]*x_start[0] + A[3][1]*x_start[2] + A[4][1]*x_start[3]
    
    j = 0
    for i in range(2, N-2):
        x[i] = A[0][i]*x_start[j] + A[1][i]*x_start[j+1]+A[3][i]*x_start[j+3]+D[i]*(-1)*b[i]
        j+=1
    
    x[N-2] = A[0][N-2]*x_start[N-4] + A[1][N-2]*x_start[N-3] + D[N-2]*(-1)*b[N-2] + A[3][N-2]*x_start[N-1]
    x[N-1] = A[0][N-1]*x_start[N-3] + A[1][N-1]*x_start[N-2] + D[N-1]*(-1)*b[N-1]
    if np.abs(result[0] - x[0]) < epsilon:
        return x
    else:
        metodaJacobiego(N, x, epsilon, result, A, D, b)


def metodaGaussaSiedlera(N, x_start, epsion, result, N_macierz, M_macierz):
    #M^-1

    d = np.zeros(N)
    x = np.zeros(N)
    d[0] = M_Macierz[0][1]
    M_odwr = np.zeros((3, N))
    for i in range(N):
        M_odwr[0][i] = 1 / M_macierz[0][i]
        for j in range(1, 3):
            M_odwr[j, i] = 0
            for k in range(i):
                M_odwr[j, i] -= M_macierz[j, k] * M_odwr[k, i]
            M_odwr[j, i] /= M_macierz[j, i]


    #M^-1*b
    x = [0.0]*N
    Mb = np.zeros(3)
    for i in range(3):
        Mb[i] = 0
        for k in range(N):
            Mb[i] += M_odwr[i, k] * b[k]

    #M^-1*N
    MN = np.zeros((3, N))
    for i in range(3):
        for j in range(N):
            if i == 0:
                MN[i, j] = M_odwr[i, j]
            else:
                MN[i, j] = M_odwr[i, j] - 0.15 * M_odwr[i - 1, j]
    j = 0

    for i in range(2, N-2):
        x[i] = MN[0, i]*x_start[j] + MN[1, i]*x_start[j+1]+MN[2, i]*x_start[j+2]+Mb[i]*b[i]
        j+=1

