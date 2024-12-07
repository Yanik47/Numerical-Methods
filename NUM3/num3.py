import time
import numpy as np
from scipy import linalg
from scipy.linalg import det
import mpmath
import matplotlib.pyplot as plt

mpmath.mp.dps = 50
#Zmienne dla sprawdzania dzialania
N = 124
razy = 80
#----------------------------------

srednia = 500
value = 0.0
wektor = [0.0] * N

def forward_substitution(juz_nie_slownik, b, n):

    for i in range(1, n):
        b[i] -= juz_nie_slownik[0][i]*b[i-1]
    return b

def backward_substitution(juz_nie_slownik, x, n):

    x[n - 1] /= juz_nie_slownik[1][n - 1]
    x[n - 2] = (x[n - 2] - juz_nie_slownik[2][n - 2] * x[n - 1]) / juz_nie_slownik[1][n - 2]


    for i in range(n - 3, -1, -1):
        x[i] = (x[i] - juz_nie_slownik[3][i] * x[i + 2] - juz_nie_slownik[2][i] * x[i + 1]) / juz_nie_slownik[1][i]
    return x

def lu(juz_nie_slownik, n):

    for i in range(1, n-2):
        juz_nie_slownik[0][i] /= juz_nie_slownik[1][i - 1]
        juz_nie_slownik[1][i] -= juz_nie_slownik[0][i] * juz_nie_slownik[2][i - 1]
        juz_nie_slownik[2][i] -= juz_nie_slownik[0][i] * juz_nie_slownik[3][i - 1]

    juz_nie_slownik[0][n-2] /= juz_nie_slownik[1][n-3]
    juz_nie_slownik[1][n-2] -= juz_nie_slownik[0][n-2] * juz_nie_slownik[2][n-3]
    juz_nie_slownik[2][n-2] -= juz_nie_slownik[0][n-2] * juz_nie_slownik[3][n-3]

    juz_nie_slownik[0][n-1] /= juz_nie_slownik[1][n-2]
    juz_nie_slownik[1][n-1] -= juz_nie_slownik[0][n-1] * juz_nie_slownik[2][n-2]

    return juz_nie_slownik
#inicialicja macierzy do sprawdzania wyniku
macierz = np.zeros((N, N))

for i in range(N):
    for j in range(0, N):
        if(i == 0):
            macierz[0][0] = 1.2
            macierz[0][1] = 0.1
            macierz[0][2] = 0.15
        else:
            if(i == j+1):
                macierz[i][j] = 0.2
            elif(i == j):
                macierz[i][j] = 1.2
            elif(i == j-1):
                macierz[i][j] = 0.1/float(i+1)
            elif(i == j-2):
                macierz[i][j] = 0.15/float((i+1)*(i+1))
    wektor[i] = float(i+1)

print("Macierz incjalizowana jako 2-wymiarowa tablica: ")
print(macierz)
print("Wektor x: ")
print(wektor)

#metoda standardowa
start = time.time()
wyznacznik_strt = det(macierz)
result = linalg.solve(macierz, wektor)
stop = time.time()
wynik = float(stop - start)

print("Czas, za ktory rozwiazany byl uklad rownan: ")
print(wynik)
print("Wyznacznik liczony przez det(A): ")
print(wyznacznik_strt)
print("Wynik(wektor y) liczony z linalg.solve: ")
print(result)

s = 0
wyznacznik = 1.0
wynik = 0.0
#slownik = {}

LU_macierz = np.zeros((N, 4))
n_lista = []
wynik_lista = []
wyz_lista = []
z_lista = []

#LU faktoryzacja
while(N<124*200):


    # = np.zeros((N, 4))
    #LU_macierz = np.zeros((N, 4))
    z = [0.0]*N
    wektor1 = [0.0]*N
    for o in range(razy):

        juz_nie_slownik = []
        #LU_macierz = np.zeros((4, N))
        #wektor1 = [0.0] * N
        #z = [0.0] * N
        wynik = 0.0
        wyznacznik = 1.0


            # Uzupełnienie diagonali macierzy A
        juz_nie_slownik.append([0] + [0.2] * (N - 1))
        juz_nie_slownik.append([1.2] * N)
        juz_nie_slownik.append([0.1 / i for i in range(1, N)] + [0])
        juz_nie_slownik.append([0.15 / i ** 2 for i in range(1, N - 1)] + [0] + [0])

        for i in range(0, N):
            wektor1[i] = i+1
        print(wektor1)


        #pomiar czasu
        start = time.time()
        LU_macierz = lu(juz_nie_slownik, N)
        y = forward_substitution(LU_macierz, wektor1, N)
        z = backward_substitution(LU_macierz, y, N)
        stop = time.time()

        wynik += stop-start



    for s in range(N):
        wyznacznik *= LU_macierz[1][s]
    wyz_lista.append(wyznacznik)
    wynik /= razy
    n_lista.append(N)
    wynik_lista.append(wynik)
    z_lista.append(z)
    N+=124
    print(wynik_lista)

print("Czas, za ktory rozwiazany zostaje uklad rownan N_min: ")
print(wynik_lista[0])
print("Wyznacznik liczony z macierzy U, N_min: ")
print(wyz_lista[0])
print("Wynik(wektor y) liczony z własnej metody: ")
print(z_lista[0])

plt.plot(np.array(n_lista), np.array(wynik_lista))
plt.xlabel("N")
plt.ylabel("Czas, potrzebny na obliczenie")
plt.show()
