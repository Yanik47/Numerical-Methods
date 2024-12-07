import numpy as np
import sympy as sp
import mpmath
import matplotlib.pyplot as plt

mpmath.mp.dps = 50

def Dh1_f(x1, h, funk):
    wynik = sp.N((funk.subs(x_symbols, x1 + h) - funk.subs(x_symbols, x1)) / h)

    return wynik

def Dh2_f(x1, h, funk):
    wynik = sp.N((funk.subs(x_symbols, x1 + h) - funk.subs(x_symbols, x1 - h)) / (2 * h))

    return wynik

#h_min = 10e-16
#h = float(0.1)
#x_value = 0.2

h_min = float(input("Prosze podac dolny prog h: \n"))

opcja_wykres = input("Prosze wybrac opcje do rysowania wykresu(wpisz jedna z 4 liczb: \n"
                     "1 - 2 wykresy typu double\n"
                     "2 - 2 wykresy typu float \n"
                     "3 - float i double dla podpunktu a\n"
                     "4 - float i double dla podpunktu b \n")
h = float(input("Prosze wprowadzic h: "))
x_value_d = float(input("Prosze wprowadzic x: "))

h_f = np.float32(h)
x_value_f = np.float32(x_value_d)

print("Wprowadzone dane: "+ str(h) + " " + str(h_min) + " " + str(x_value_d) + "\n")

wybor_f = input("Prosze wybrac funkcje do dzialania: \n")

x_symbols = sp.symbols('x')
funk = sp.sympify(wybor_f)
poch = funk.diff(x_symbols)
poch_result_d = poch.subs(x_symbols, x_value_d)
poch_result_f = poch.subs(x_symbols, x_value_f)

blad1 = []
blad2 = []
h_lista = []

while(h > h_min):

    h_d = h
    h_f = np.float32(h)
    y1 = Dh1_f(x_value_d, h_d, funk)
    print("wynik1(double): " + str(y1))
    y2 = Dh1_f(x_value_f, h_f, funk)
    print("wynik1(float): " + str(y2))
    y3 = Dh2_f(x_value_d, h_d, funk)
    print("wynik2(double): " + str(y3))
    y4 = Dh2_f(x_value_f, h_f, funk)
    print("wynik2(float): " + str(y4))

    if (opcja_wykres == "2"):
        blad1.append(np.float32(abs(poch_result_f - y2)))
        blad2.append(np.float32(abs(poch_result_f - y4)))
    elif(opcja_wykres == "1"):
        blad1.append(float(abs(poch_result_d - y1)))
        blad2.append(float(abs(poch_result_d - y3)))
    elif(opcja_wykres == "3"):
        blad1.append(np.float32(abs(poch_result_f - y1)))
        blad2.append(float(abs(poch_result_d - y2)))
    elif(opcja_wykres == "4"):
        blad1.append(np.float32(abs(poch_result_f - y4)))
        blad2.append(float(abs(poch_result_d - y3)))
    else:
        print("Wprowadzono nieznana opcje")

    result_1 = blad1[len(blad1)-1]
    result_2 = blad2[len(blad2)-1]

    h = h/1.1
    h_lista.append(h)

    print("Znaczenie pochodnej nie ze wzoru: " +str(poch_result_d)+ "\n")
    print("Znaczenie pochodnej pierwszy wzor(double): " + str(y1) + "\n")
    print("Znaczenie pochodnej pierwszy wzor(float): " + str(y2) + "\n")
    print("          ")
    print("Blad obliczen pierwszy wzor: " + str(result_1))
    print("\n")

    print("Znaczenie pochodnej drugi wzor(double): " + str(y3)+ "\n")
    print("Znaczenie pochodnej pierwszy wzor(float): " + str(y4) + "\n")
    print("          ")
    print("Blad obliczen drugi wzor:"+ str(result_2))
    print("\n")

plt.plot(np.log10(np.array(h_lista)), np.log10(np.array(blad1)))
plt.plot(np.log10(np.array(h_lista)), np.log10(np.array(blad2)))
plt.xlabel("h")
plt.ylabel("Blad przyblizenia")
plt.show()