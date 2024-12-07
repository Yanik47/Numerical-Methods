import time
import numpy as np
from scipy import linalg
from scipy.linalg import det
import mpmath
import matplotlib.pyplot as plt

mpmath.mp.dps = 50
#Zmienne dla sprawdzania dzialania
N = 124
razy = 30
#----------------------------------

srednia = 500
value = 0.0
wektor = [0.0] * N

def forward_substitutionn(slownik, b):
    n = len(b)
    x = [0]*n

    for i in range(n):
        x[i] = b[i]
        for j in range(i):
            x[i] -= get_value(slownik, i, j) * x[j]
    return x

def backward_substitution(tab, x):
    n = len(x)
    y = [0.0]*n
    for i in range( n-1, -1, -1):
        y[i] = x[i]
        for j in range(i+1, n):
            y[i] -= x[j]*get_value(tab, i, j)
        x[i] = y[i]/get_value(tab, i, i)
    return x

def lu(slownik, n):
    set_value(slownik, 1, 0, get_value(slownik, 1, 0) / get_value(slownik, 0, 0))
    for i in range(1, n):
        set_value(slownik, i, i, get_value(slownik, i, i) - get_value(slownik, i, i - 1) * get_value(slownik, i - 1, i))
        set_value(slownik, i + 1, i, get_value(slownik, i + 1, i) / get_value(slownik, i, i))
        set_value(slownik, i, i + 1,
                  get_value(slownik, i, i + 1) - get_value(slownik, i, i - 1) * get_value(slownik, i - 1, i + 1))
    return slownik

def set_value(tab, row, col, value):
    tab[(row, col)] = value

def get_value(tab, row, col):
    return tab.get((row, col), 0)

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
slownik = {}
LU_macierz = {}

n_lista = []
wynik_lista = []
wyz_lista = []
z_lista = []

#LU faktoryzacja
while(N<1488*2):
    wektor1 = [0.0]*N
    z = [0.0] * N
    for o in range(razy):
        slownik.clear()
        wynik = 0.0
        wyznacznik = 1.0

        for i in range(0, N):
            #j = 0
            for j in range(0, N):
                if (i == j + 1):
                    value = 0.2
                    set_value(slownik, i, j, value)
                elif (i == j):
                    value = 1.2
                    set_value(slownik, i, j, value)
                elif (i == j - 1):
                    value = 0.1 / float(i + 1)
                    set_value(slownik, i, j, value)
                elif (i == j - 2):
                    value = 0.15 / float((i + 1) * (i + 1))
                    set_value(slownik, i, j, value)
            wektor1[i] = i+1


        #pomiar czasu
        start = time.time()
        LU_macierz = lu(slownik, N)
        y = forward_substitutionn(LU_macierz, wektor1)
        z = backward_substitution(LU_macierz, y)
        stop = time.time()

        wynik += stop-start
    for s in range(N):
        wyznacznik *= get_value(LU_macierz, s, s)
    wyz_lista.append(wyznacznik)
    wynik /= razy
    n_lista.append(N)
    wynik_lista.append(wynik)
    z_lista.append(z)
    N+=248
    print(wynik_lista)

print("Czas, za ktory rozwiazany zostaje uklad rownan N_min: ")
print(wynik_lista[0])
print("Wyznacznik liczony z macierzy U, N_min: ")
print(wyz_lista[0])
print("Wynik(wektor y) liczony z własnej metody: ")
print(z_lista[0])

plt.plot(np.array(n_lista), np.array(wynik_lista))
plt.xlabel("(log10) N")
plt.ylabel("(log10) Czas, potrzebny na obliczenie")
plt.show()


