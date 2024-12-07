
import time
import numpy as np
import matplotlib.pyplot as plt
N = 80
np.set_printoptions(precision=15, suppress=True)
time_list = []
x_list = []
N_list = []
time_b_list = []
result_b_list = []
def moja_piekna_funckja(N, razy, rozmiar):
    while(N<=80*rozmiar):
        r = 0
        wynik = 0.0
        #####################
        wynik_b = 0.0
        result = 0.0
        #####################
        #while(r<razy):
        A_b = np.zeros((N, N))
        wektor_b = [5.0] * N
        #--------------BIBLIOTECZNA METODA--------------------------------------------
        for i in range(N):
            for j in range(N):
                if (i == j):
                    A_b[i][j] = 12
                elif (j == i + 1):
                    A_b[i][j] = 8
                else:
                    A_b[i][j] = 1
        start_b = time.time()
        result = np.linalg.solve(A_b, wektor_b)
        stop_b = time.time() - start_b
        wynik_b +=stop_b
        while(r<razy):
        #--------------METODA WLASNA--------------------------------------------------
            x = y = z = np.zeros(N)
            Ann = 11.0
            Ann1 = 7.0
            b = 5.0
            start = time.time()
            z[len(z)-1] = b/Ann
            y[len(z)-1] = 1/Ann
            for i in range(N-2, -1, -1):
                z[i] = (b - z[i+1]*Ann1)/Ann
                y[i] = (1 - y[i+1]*Ann1)/Ann

            skalar = sum(z)/(1+sum(y))

            for i in range(N):
                x[i] = (z[i] - skalar*y[i])*b
                #print(x)
            stop = time.time() - start
            r+=1
            wynik+=stop
        #----------------------------------------------------------------------------

        wynik /=razy
        #wynik_b /= razy

        ##################################
        result_b_list.append(result)
        time_b_list.append(wynik_b)
        ##################################
        N_list.append(N)
        N += 150
        time_list.append(wynik)
        x_list.append(x)

    print("Wynik metoda wlasna: \n")
    print(x_list[0])
    print("Wynik metoda biblioteczna: \n")
    print(result_b_list[0])
    print("Wynik pomiaru czasu z wlasnej metody: \n")
    print(time_list[0])
    print("Wynik pomiaru czasu z metudy bibliotecznej: \n")
    print(time_b_list[0])
    plt.plot(np.array(N_list), np.array(time_list))
    #########################################################
    #plt.plot(np.array(N_list), np.array(time_b_list))
    #########################################################

    plt.xlabel("N")
    plt.ylabel("Czas, potrzebny na obliczenie")
    plt.show()
moja_piekna_funckja(80, 20, 50)
#moja_piekna_funckja(80, 40, 50)
#moja_piekna_funckja(80, 60, 50)
#moja_piekna_funckja(80, 80, 50)
