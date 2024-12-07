import numpy as np
import matplotlib.pyplot as plt


def wektor_wlasny(eigenvalue, A, max_iterations=1000, epsilon=1e-12):
    n = A.shape[0]
    v = np.random.rand(n)
    for _ in range(max_iterations):
        v_new = np.linalg.solve(A -eigenvalue *np.eye(n), v)
        v_new = v_new / np.linalg.norm(v_new)
        if np.linalg.norm(v_new - v) < epsilon:
            break
        v = v_new
    return v

def wilkinson_shift(A):
    n = A.shape[0]
    alpha = A[n-1, n-1]
    beta = A[n-1, n-2]

    B = np.array([[alpha, beta], [beta, alpha]])
    mu, _ = np.linalg.eig(B)

    shift = mu[0] if np.abs(mu[0] - alpha) < np.abs(mu[1] - alpha) else mu[1]

    return shift

M = np.array([[8.0, 1.0, 0.0, 0.0],
              [1.0, 7.0, 2.0, 0.0],
              [0.0, 2.0, 6.0, 3.0],
              [0.0, 0.0, 3.0, 5.0]])

#biblioteczna metoda
values, vectors = np.linalg.eig(M)
index_max_value = np.argmax(np.abs(values))
max_value = values[index_max_value]
max_value_vector = vectors[:, index_max_value]

#potegowa metoda
w = [1, 1, 1, 1]
z = [1, 1, 1, 1]
e = 1.0
lam = 0.0
lam_pot_lista = []
iter_pot_lista = []
epsilon = 10e-13
i = 0
I = np.eye(4)
#prosze zakomentowac dla generowania wykresu bez przesuniecia
M[0, 0] = M[0, 0] - (2.07816482+8.14771771)/2.0
M[1, 1] = M[1, 1] - (2.07816482+8.14771771)/2.0
M[2, 2] = M[2, 2] - (2.07816482+8.14771771)/2.0
M[3, 3] = M[3, 3] - (2.07816482+8.14771771)/2.0
##############################################################
while e > epsilon and i < 5000:
    w[0] = M[0, 0]*z[0] + M[0, 1]*z[1]
    w[1] = M[1, 0]*z[0] + M[1, 1]*z[1] + M[1, 2]*z[2]
    w[2] = M[2, 1]*z[1] + M[2, 2]*z[2] + M[2, 3]*z[3]
    w[3] = M[3, 2]*z[2] + M[3, 3]*z[3]
    lam = np.max(np.abs(w))
    e = np.sqrt((w[0] - lam*z[0])**2 + (w[1] - lam*z[1])**2 + (w[2] - lam*z[2])**2 + (w[3]-lam*z[3])**2)
    z = w/lam
    #prosze zakomentowac dla generowania wykresu bez przesuniecia
    ############################
    lam += (2.07816482+8.14771771)/2.0
    ###########################
    lam_pot_lista.append(np.abs(max_value - lam))
    i += 1
    iter_pot_lista.append(i)
    #z - wektor wlasny
M = np.array([[8.0, 1.0, 0.0, 0.0],
              [1.0, 7.0, 2.0, 0.0],
              [0.0, 2.0, 6.0, 3.0],
              [0.0, 0.0, 3.0, 5.0]])
v0 = wektor_wlasny(lam, M)

#metoda QR
QT = np.zeros((4, 4))
Q = np.zeros((4, 4))
R = np.zeros((4, 4))
M_iter = np.zeros((4, 4))
PT = np.eye(4)
P = np.eye(4)

QR_wykres_1 = []
QR_wykres_2 = []
QR_wykres_3 = []
QR_wykres_4 = []
iter_list = []
j = 0

#M_low =np.array([[6, 3],
#                [3, 5]])
#values_low, vectors_low = np.linalg.eig(M_low)
#M = np.add(M, -P*values_low[1])

while np.abs(M[2, 1]) > epsilon and np.abs(M[1, 0]) > epsilon and np.abs(M[3, 2]) > epsilon and j < 150:
    Q, R = np.linalg.qr(M_iter)
    QT = np.transpose(Q)
    PT = np.dot(QT, PT)
    P = np.dot(P, Q)
    M_iter = np.dot(PT, np.dot(M, P))

    QR_wykres_1.append(np.abs(M_iter[0, 0] - values[0]))
    QR_wykres_2.append(np.abs(M_iter[1, 1] - values[1]))
    QR_wykres_3.append(np.abs(M_iter[2, 2] - values[2]))
    QR_wykres_4.append(np.abs(M_iter[3, 3] - values[3]))

    j += 1
    iter_list.append(j)
v = np.zeros((4, 4))

print("Macierz jest postaci: \n")
print(M_iter)

#wyniki
for i in range(4):
    v[i] = wektor_wlasny(M_iter[i, i], M)

    print("To jest wartosc wlasna macierzy M obliczony metoda QR: \n")
    print(M_iter[i, i])
    print("\n To jest wektor wlasny macierzy M obliczony metoda QR: \n")
    print(v[i])

print("To jest najwieksza wartosc wlasna obliczona metoda iteracyjna: \n")
print(lam)
print("Odpowiedni wektor: ")
print(v0)

print("To jest najwieksza wartosc wlasna obliczona metoda biblioteczna: \n")
print(max_value)
print("Odpowiedni wektor: \n")
print(max_value_vector)

#tworzenie wykresow dla potegowej metody
plt.xlabel('Iteracja(i)')
plt.grid(True)
plt.ylabel("log10|eigenValue(iteracyjna) - eigenValue(dokladna)|")
plt.plot(np.array(iter_pot_lista), np.log10(np.array(lam_pot_lista)))
plt.title('Najwieksza wartosc wlasna dla kolejnej iteracji metody potegowej')
plt.show()

#tworzenie wykresow dla metody QR
plt.xlabel('Iteracja(i)')
plt.grid(True)
plt.ylabel("log10|eigenValue(iteracyjna) - eigenValue(dokladna)|")
plt.plot((np.array(iter_list)), np.log10(np.array(QR_wykres_1)))
plt.plot((np.array(iter_list)), np.log10(np.array(QR_wykres_2)))
plt.plot((np.array(iter_list)), np.log10(np.array(QR_wykres_3)))
plt.plot(((np.array(iter_list))), np.log10(np.array(QR_wykres_4)))
plt.legend(['Pierwsza wartosc wlasna', 'Druga wartosc wlasna', 'Trzecia wartosc wlasna', 'Czwarta wartosc wlasna'])
plt.title('Wartosci wlasne dla kolejnych iteracji metody QR')
plt.show()
