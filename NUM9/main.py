import math
import matplotlib.pyplot as plt
import numpy as np

epsilon = 10e-8
x11 = 0
x22 = math.pi/2
iteracja = 0
max_iter = 50
my_function = {1: "arcsin(x)-0.4", 2: "(arcsin(x)-0.4)^2"}
err = []


def bisekcji(x1, x2, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter and f(x1, ftype) * f(x2, ftype) < 0:
        iteracja += 1
        x3 = (x1 + x2) / 2
            #print(x3)
        err.append(abs(x3 - math.asin(0.4)))
        if abs(x3 - x2) <= epsilon:
                #print("hi_ bis")
            return x3
        else:
            if f(x1, ftype) * f(x3, ftype) < 0:
                return bisekcji(x1, x3, ftype)
            elif f(x3, ftype) * f(x2, ftype) < 0:
                return bisekcji(x3, x2, ftype)
            else:
                return x3
    else:
        return


def falsi(x1, x2, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter and f(x1, ftype) * f(x2, ftype) < 0:
        iteracja += 1
        x3 = (f(x1, ftype) * x2 - f(x2, ftype) * x1) / (f(x1, ftype) - f(x2, ftype))
            #print(x3)
        err.append(abs(x3 - math.asin(0.4)))
        if abs(x3 - x2) <= epsilon:
                #print("hi")
            return x3
        else:
            if f(x1, ftype) * f(x3, ftype) < 0:
                return falsi(x1, x3, ftype)
            elif f(x3, ftype) * f(x2, ftype) < 0:
                return falsi(x3, x2, ftype)
            else:
                return x3
    else:
        return


def siecznych(x1, x2, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter:
        iteracja += 1
        x3 = (f(x1, ftype) * x2 - f(x2, ftype) * x1) / (f(x1, ftype) - f(x2, ftype))
        err.append(abs(x3 - math.asin(0.4)))
        #print("x3: " + str(x3) + " x2: " + str(x1))
        if abs(x3 - x1) < epsilon:
            return x3
        else:
            return siecznych(x3, x1, ftype)
    else:
        return

def g_siecznych(x1, x2, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter:
        iteracja += 1
            #x3 = (f(x1, ftype) * x2 - f(x2, ftype) * x1) / (f(x1, ftype) - f(x2, ftype))
        ux1 = f(x1, ftype)/(2 * math.cos(x1)*(math.sin(x1) - 2/5))
        ux2 = f(x2, ftype)/(2 * math.cos(x2)*(math.sin(x2) - 2/5))
        x3 = (ux1*x2 - ux2*x1)/(ux1 - ux2)
        #print("x3: " + str(x3))
        #print("x2: " + str(x2))
        #print("x1: " + str(x1))
        err.append(abs(x3 - math.asin(0.4)))
            # print("x3: " + str(x3) + " x2: " + str(x1))
        if abs(x3 - x2) < epsilon:
            return x3
        else:
            return g_siecznych(x3, x1, ftype)
    else:
        return


def f_newtona(x, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter:
        iteracja += 1
        df = forward_difference(f, x, ftype)
        #df = math.cos(x)
        ux = f(x, ftype)/df
        x1 = x - ux
        #print(x1)
        err.append(abs(x1 - math.asin(0.4)))
        if abs(x1 - x) <= epsilon:
            return x1
        else:
            return f_newtona(x1, ftype)
    return


#??????????????????
def g_newtona(x, ftype):
    global iteracja
    global max_iter
    if iteracja < max_iter:
        iteracja +=1
        # df1 = forward_difference(f, x, ftype)
        df1 = 2 * math.cos(x)*(math.sin(x) - 2/5)
        #df2 = forward_difference(df1, x, ftype)
        df2 = 2 * (math.cos(x)**2 - (math.sin(x)-2/5)*math.sin(x))
        print("x: " + str(x))
        u1x = f(x, ftype)/df1
        u2x = 1 - ((f(x, ftype)/df1**2)*df2)
        x1 = x - u1x/u2x
        #print(x1)
        err.append(abs(x1 - math.asin(0.4)))
        if(abs(x1-x) <= epsilon):
            return x1
        else:
            return g_newtona(x1, ftype)


def f(x, ftype):
    if ftype == 1:
        return math.sin(x) - 0.4
    else:
        return (math.sin(x) - 0.4) ** 2



def forward_difference(f, x, ftype):
    h = 10e-7
    df = (f(x + h, ftype) - f(x, ftype)) / h
    return df



#print("fds")
result_bi = bisekcji(x11, x22, 1)
#print("fdsf")
print("Wynik metody bisekcji: " + str(result_bi))
if result_bi is None:
    print("Nie udalo sie obliczyc metoda bisekcji dla funkcji: ")
    print(my_function.get(1))
else:
    result_iter_bi = iteracja
    print("Ilosc iteracji dla metody bisekcji: " + str(result_iter_bi))
    iter_array = [i for i in range(1, result_iter_bi+1)]
    plt.plot(iter_array, np.log10(np.array(err)), marker='o')
#####################################
#plt.title('Wykresy bledu obliczeniowego dla kolejnej iteracji funkcji ' + my_function.get(1))
#plt.xlabel('Iteracja')
#plt.ylabel('x(i) - x(dokladny)')
#plt.grid()
#plt.legend(['Metoda bisekcji'])
#plt.savefig("wykres3.svg")
#plt.show()
#plt.clf()
########################################
iteracja = 0
err = []

result_fa = falsi(x11, x22, 1)
print("Wynik metody falsi: " + str(result_fa))
#print(iteracja)
if result_fa is None:
    print("Nie udalo sie obliczyc metoda falsi dla funkcji:  ")
    print(my_function.get(1))
else:
    result_iter_fa = iteracja
    print("Ilosc iteracji dla metody falsi: " + str(result_iter_fa))
    iter_array = [i for i in range(1, result_iter_fa+1)]
    plt.plot(iter_array, np.log10(np.array(err)), marker='o')
    # plot
#############################################
#plt.title('Wykresy bledu obliczeniowego dla kolejnej iteracji funkcji ' + my_function.get(1))
#plt.xlabel('Iteracja')
#plt.ylabel('x(i) - x(dokladny)')
#plt.grid()
#plt.legend(['Metoda regula falsi'])
#plt.savefig("wykres4.svg")
#plt.show()
#plt.clf()
############################################
iteracja = 0
err = []

result_sie = siecznych(x11, x22, 1)
print("Wynik metody siecznych: " + str(result_sie))
#print(iteracja)
if result_sie is None:
    print("Nie udalo sie obliczyc metoda siecznych dla funkcji: ")
    print(my_function.get(1))
else:
    result_iter_sie = iteracja
    print("Ilosc iteracji dla metody siecznych: " + str(result_iter_sie))
    iter_array = [i for i in range(1, result_iter_sie+1)]
    plt.plot(iter_array, np.log10(np.array(err)),marker='o')
    # plot
##############################################################
#plt.title('Wykresy bledu obliczeniowego dla kolejnej iteracji funkcji ' + my_function.get(1))
#plt.xlabel('Iteracja')
#plt.ylabel('x(i) - x(dokladny)')
#plt.grid()
#plt.legend(['Metoda siecznych'])
#plt.savefig("wykres5.svg")
#plt.show()
#plt.clf()
############################################################33
iteracja = 0
err = []

result_New = f_newtona(x11, 1)
print("Wynik metody Newtona: " + str(result_New))
if result_New is None:
    print("Nie udalo sie obliczyc metoda Newtona dla funkcji: ")
    print(my_function.get(1))
else:
    result_iter_New = iteracja
    print("Ilosc iteracji dla metody Newtona: " + str(result_iter_New))
    iter_array = [i for i in range(1, result_iter_New+1)]
    plt.plot(iter_array, np.log10(np.array(err)), marker='o')
    # plot - iteracja na x, err na y

iteracja = 0
err = []
##############################################################################
plt.title('Wykresy bledu obliczeniowego dla kolejnej iteracji funkcji ' + my_function.get(1))
plt.xlabel('Iteracja')
plt.ylabel('x(i) - x(dokladny)')
plt.grid()
plt.legend(['Metoda Bisekcji','Metoda Falsi','Metoda Siecznych','Metoda Newtona'])
plt.savefig("wykres6.svg")
plt.show()
plt.clf()
#############################################################################
result_sie = g_siecznych(x11, x22, 2)
print("Wynik metody siecznych: " + str(result_sie))
#print(iteracja)
if result_sie is None:
    print("Nie udalo sie obliczyc metoda siecznych dla funkcji: ")
    print(my_function.get(2))
else:
    result_iter_sie = iteracja
    print("Ilosc iteracji dla metody siecznych: " + str(result_iter_sie))
    iter_array = [i for i in range(1, result_iter_sie+1)]
    plt.plot(iter_array, np.log10(np.array(err)), marker='o')
    # plot

iteracja = 0
err = []

result_New = g_newtona(x11, 2)
print("Wynik metody Newtona: " + str(result_New))
if result_New is None:
    print("Nie udalo sie obliczyc metoda Newtona dla funkcji: ")
    print(my_function.get(2))
else:
    result_iter_New = iteracja
    print("Ilosc iteracji dla metody Newtona: " + str(result_iter_New))
    iter_array = [i for i in range(1, result_iter_New+1)]
    plt.plot(iter_array, np.log10(np.array(err)), marker='o')
    # plot - iteracja na x, err na y

plt.title('Wykresy bledu obliczeniowego dla kolejnej iteracji funkcji ' + my_function.get(2))
plt.xlabel('Iteracja')
plt.ylabel('x(i) - x(dokladny)')
plt.grid()
plt.legend(['Metoda siecznych', 'Metoda Newtona'])
plt.savefig("wykres2.svg")
plt.show()
plt.clf()