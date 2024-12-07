import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

x = []
x_ale_zmienna = -1.0
##################wartosci na zmiane przy potrzebie recznej zmiany########
odleglosc_punktow = 0.01
##########################################################################
while x_ale_zmienna <= 1.0:
    x.append(x_ale_zmienna)
    x_ale_zmienna += odleglosc_punktow


def w_calculation(y, x, n, len_x, x_i):
    w = []
    for i in range(len_x + 1):
        temp = 0
        for k in range(n):
            phi = 1
            for j in range(n):
                if j != k:
                    phi *= (x[i] - x_i[j]) / (x_i[k] - x_i[j])
            temp += y[k] * phi
        w.append(temp)
    return w


def wzor_a(n):
    x_a = []
    for i in range(n):
        x_a.append(-1.0 + 2.0 * (i / (n - 1.0)))
    print(x_a)
    return x_a


def wzor_b(n):
    x_b = []
    for i in range(n):
        x_b.append(np.cos((((2.0 * i) + 1.0) / (2.0 * (n + 1.0))) * np.pi))
    return x_b


def first_function(n, x_i):
    y = []
    for i in range(n):
        y.append(1.0 / (1.0 + 50.0 * (x_i[i] * x_i[i])))
    print(y)
    return y


def second_function(n, x_i):
    y = []
    for i in range(n):
        y.append(np.cos(5 * np.pi * x_i[i] + np.sin(3 * np.pi * x_i[i])))
    return y


def third_function(n, x_i):
    y = []
    for i in range(n):
        y.append(5 * pow(x_i[i], 5) - 4 * x_i[i] ** 5 + 3 * x_i[i] ** 3 - 2 * x_i[i] ** 2 + x_i[i])
    return y


def bib_interpolation(x_i, y, x):
    x_i = np.array(x_i)
    y = np.array(y)
    min_x, max_x = min(x_i), max(x_i)
    valid_indices = (x_i >= min_x) & (x_i <= max_x)
    x_i_valid = x_i[valid_indices]
    interpolation_function_a = interp1d(x_i_valid, y[valid_indices], kind='cubic')
    x_interp = np.linspace(min_x, max_x, 1000)
    y_interp = interpolation_function_a(x_interp)

    # Rysuj wykres
    plt.plot(x_interp, y_interp, linestyle='dashed', linewidth='5')


def plot1(x, n, wzor, opcja):
    x_i = wzor
    y_norm = []
    if (opcja == 1):
        y = first_function(n, x_i)
        y_norm = first_function(len(x), x)
    elif (opcja == 2):
        y = second_function(n, x_i)
        y_norm = second_function(len(x), x)
    elif (opcja == 3):
        y = third_function(n, x_i)
        y_norm = third_function(len(x), x)
    else:
        return

    w = w_calculation(y, x, n, len(x) - 1, x_i)
    y_vals = []
    for i in range(len(x) - 1):
        y_vals.append(np.abs(y_norm[i] - w[i]))
    plt.plot(x, w, '-')
    return np.max(y_vals)


err = []
num_err = [1, 2, 3, 4, 5]

# pierwsza funckja a
plt.xlim(-1, 1)
plt.ylim(-5, 5)
bib_interpolation(wzor_a(20), first_function(20, wzor_a(20)), x)
err.append(plot1(x, 3, wzor_a(3), 1))
err.append(plot1(x, 4, wzor_a(4), 1))
err.append(plot1(x, 5, wzor_a(5), 1))
err.append(plot1(x, 7, wzor_a(7), 1))
err.append(plot1(x, 30, wzor_a(30), 1))

plt.title('Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=\frac{1}{1+50x^2}$ i siatki $x_i=-1 + 2\frac{i}{n+1}$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_3(x)', 'W_4(x)', 'W_5(x)', 'W_7(x)', 'W_30(x)'])
plt.savefig("wykres.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 1000)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err.svg")
plt.show()

plt.xlim(-1, 1)
plt.ylim(-5, 5)
err = []
bib_interpolation(wzor_b(20), first_function(20, wzor_b(20)), x)
err.append(plot1(x, 3, wzor_b(3), 1))
err.append(plot1(x, 4, wzor_b(4), 1))
err.append(plot1(x, 5, wzor_b(5), 1))
err.append(plot1(x, 7, wzor_b(7), 1))
err.append(plot1(x, 30, wzor_b(30), 1))

plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_3(x)', 'W_4(x)', 'W_5(x)', 'W_7(x)', 'W_30(x)'])
plt.title('Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=\frac{1}{1+50x^2}$ i siatki $x_i=\cos(\pi\frac{2i+1}{2('
          r'n+1)})$')
plt.savefig("wykres1.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 10)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err1.svg")

plt.show()
plt.clf()

plt.xlim(-1, 1)
plt.ylim(-5, 5)
err = []
bib_interpolation(wzor_a(20), second_function(20, wzor_a(20)), x)
err.append(plot1(x, 3, wzor_a(3), 2))
err.append(plot1(x, 4, wzor_a(4), 2))
err.append(plot1(x, 6, wzor_a(6), 2))
err.append(plot1(x, 7, wzor_a(7), 2))
err.append(plot1(x, 25, wzor_a(25), 2))

plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_3(x)', 'W_4(x)', 'W_6(x)', 'W_7(x)', 'W_25(x)'])
plt.title(
    'Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=\cos(\pi5x) + sin(\pi3x)$ i siatki $x_i=-1 + 2\frac{i}{n+1}$')
plt.savefig("wykres2.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 1000)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err2.svg")

plt.show()
plt.clf()

plt.xlim(-1, 1)
plt.ylim(-5, 5)
err = []
bib_interpolation(wzor_b(20), second_function(20, wzor_b(20)), x)
err.append(plot1(x, 3, wzor_b(3), 2))
err.append(plot1(x, 4, wzor_b(4), 2))
err.append(plot1(x, 6, wzor_b(6), 2))
err.append(plot1(x, 7, wzor_b(7), 2))
err.append(plot1(x, 25, wzor_b(25), 2))
plt.title('Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=\cos(\pi5x) + sin(\pi3x)$ i siatki $x_i=\cos(\pi\frac{'
          r'2i+1}{2(n+1)})$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_3(x)', 'W_4(x)', 'W_6(x)', 'W_7(x)', 'W_25(x)'])
plt.savefig("wykres3.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 20)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err3.svg")

plt.show()
plt.clf()

plt.xlim(-1, 1)
plt.ylim(-5, 5)
err = []
bib_interpolation(wzor_a(20), third_function(20, wzor_a(20)), x)
err.append(plot1(x, 2, wzor_a(2), 3))
err.append(plot1(x, 3, wzor_a(3), 3))
err.append(plot1(x, 4, wzor_a(4), 3))
err.append(plot1(x, 8, wzor_a(8), 3))
err.append(plot1(x, 15, wzor_a(15), 3))
plt.title(
    'Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=5x^5 - 4x^4 + 3x^3 - 2x^2 + x$ i siatki $x_i=-1 + 2\frac{i}{n+1}$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_2(x)', 'W_3(x)', 'W_4(x)', 'W_8(x)', 'W_15(x)'])
plt.savefig("wykres4.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 20)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err4.svg")

plt.show()
plt.clf()

plt.xlim(-1, 1)
plt.ylim(-5, 5)
err = []
bib_interpolation(wzor_b(20), third_function(20, wzor_b(20)), x)
err.append(plot1(x, 2, wzor_b(2), 3))
err.append(plot1(x, 3, wzor_b(3), 3))
err.append(plot1(x, 4, wzor_b(4), 3))
err.append(plot1(x, 8, wzor_b(8), 3))
err.append(plot1(x, 15, wzor_b(15), 3))
plt.title(
    'Wielomiany interpolacyjne dla funkcji\n' r'$f(x)=5x^5 -4x^4 + 3x^3 - 2x^2 + x$ i siatki $x_i=\cos(\pi\frac{2i+1}{2(n+1)})$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend(['Dokladny wykres z biblioteki', 'W_2(x)', 'W_3(x)', 'W_4(x)', 'W_8(x)', 'W_15(x)'])
plt.savefig("wykres5.svg")
plt.show()
plt.clf()

plt.xticks(np.arange(1, 6, 1))
plt.ylim(0, 10)
# Chebyshev
plt.title("Blad interpolacji")
plt.xlabel('Numer wykresu z legendy wyzej')
plt.ylabel('Norma Chebysheva dla tego wykresu')
plt.grid(axis='y')
plt.bar(num_err, err)
plt.savefig("err5.svg")

plt.show()
plt.clf()
