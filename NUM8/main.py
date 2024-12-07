import numpy as np
import csv
import matplotlib.pyplot as plt
#liczenie wartosci y dla funkcji
def poly_1(coefficients, x):
    return coefficients[0]*x**2 + coefficients[1]*np.sin(x) + coefficients[2]*np.cos(5*x) + coefficients[3]*np.exp(-x)

#liczenie wartosci y dla funkcji
def poly_2(cof, x):
    return cof[0]*np.sin(3*x) + cof[1]*np.cos(2*x) + cof[2]*x**2 + cof[3]*np.exp(-0.2*x)

#liczenie macierzy VanderMonde'a
def macierz_VanderMonde_1(x_data, degree):
    A = np.zeros((len(x_data), degree + 1))
    for indeks in range(len(x_data)):
        A[indeks][0] = x_data[indeks] ** 2
        A[indeks][1] = np.sin(x_data[indeks])
        A[indeks][2] = np.cos(5*x_data[indeks])
        A[indeks][3] = np.exp(-x_data[indeks])
    return A
#liczenie macierzy VanderMonde'a
def macierz_VanderMonde_2(x_data, degree):
    A = np.zeros((len(x_data), degree + 1))
    for indeks in range(len(x_data)):
        A[indeks][2] = x_data[indeks] ** 2
        A[indeks][0] = np.sin(3*x_data[indeks])
        A[indeks][1] = np.cos(2*x_data[indeks])
        A[indeks][3] = np.exp(-0.2*x_data[indeks])
    return A

def aproksymacja(x_data, y_data, degree, A):
    #rozklad SVD
    U, S, V = np.linalg.svd(A)
    # U^T
    U = U.T
    # S^+
    #liczenie macierzy pseudoodwrotnej S
    S_matrix = np.zeros((degree + 1, len(x_data)))
    for i in range(degree + 1):
        for j in range(degree + 1):
            if (j == i):
                S_matrix[i, j] = 1 / S[i]
    #liczenie macierzy pseudoodwrotnej A
    A_pseudo = np.dot(V.T, np.dot(S_matrix, U))

    coefficients = np.dot(A_pseudo, y_data)
    print("a: ", coefficients)

    return coefficients

file_path = 'dane.txt'
x_data = []
y_data = []
#wczytujemy dane z pliku
with open(file_path, 'r') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        if row:
            x_data.append(np.float64(row[0]))
            y_data.append(np.float64(row[1]))
degree = 3

A = macierz_VanderMonde_1(x_data, degree)

#metoda wlasna
coefficients = aproksymacja(x_data, y_data, degree, A)
result = []
for i in range(len(x_data)):
    result.append(poly_1(coefficients, x_data[i]))
#metoda biblioteczna
coef_bib = np.linalg.lstsq(A, y_data, rcond=None)[0]
result_bib = []
for i in range(len(x_data)):
    result_bib.append(poly_1(coef_bib, x_data[i]))


plt.plot(x_data, y_data, 'o')
plt.plot(x_data, result_bib, marker='o')
plt.plot(x_data, result, marker='*', color='r')

plt.title('Wykres dla funkcji\n' r'$F(x)=ax^2 + b\sin(x)+c\cos(5x) + \frac{d}{e^x}$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Wykres z biblioteki', 'Wlasna metoda', 'punkty podane w zadaniu'])
plt.savefig("wykres.svg")
plt.show()
plt.clf()
###PARAMETRY DO ZMIANY
a = [2, 1.5, 0.5, 0.8]
il_pkt = 100 #ilosc punktow
std_dev = 1 #zaburzenie
######################
#generowanie x oraz y
x_values = np.linspace(-5, 4, il_pkt)
y_values = poly_2(a, x_values)
#Szum
random_noise = np.random.normal(loc=0, scale=std_dev, size=il_pkt)
y_values_noice = y_values + random_noise

#Algorytm
A = macierz_VanderMonde_2(x_values, degree)
cof = aproksymacja(x_values, y_values_noice, degree, A)
result = []
for i in range(len(x_values)):
    result.append(poly_2(cof, x_values[i]))
#wykresy
plt.plot(x_values, y_values, marker='o')
plt.plot(x_values, result, marker='*')
plt.plot(x_values, y_values_noice, 'o')

plt.title('Wykres dla funkcji\n' r'$F(x)=a*sin(3x) + b\cos(2x) + cx^2) + \frac{d}{0.8e^{-0.2x}}$')
plt.xlabel('Zaburzenie jest rowne 1,Punktow jest 100                                     x')
plt.ylabel('y')
plt.grid()
plt.legend(['Wykres z biblioteki', 'Wlasna metoda', 'punkty podane z szumem'])
plt.savefig("wykres6.svg")
plt.show()
plt.clf()



