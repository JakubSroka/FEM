import sys
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

"""
zad. 9
-u'' + k(x)u' = f(x)
u(0) = u0
u'(2) = g

k(x) = 1 ,x <- [0;1)
k(x) = 2 ,x <- [1,2]
"""

def k(x):
    if x < 1:
        return 1
    else:
        return 2

def MES(u0, g, n, f):
    "obliczamy odległość między punktami"
    h = 2.0/n
    "Obliczam wyrazy macierzy głównej"
    P1 = lambda x: -1/h + 1/2 * k(x)
    P2 = 2/h
    P3 = lambda x: -1/h - 1/2 * k(x)

    "Wstepne uzupelnienie macierzy glownej poszerzonej o warunek neumana i dirichleta"
    macierz = np.zeros((n, n))
    for i in range(n):
        macierz[i, i] = P2
        if (i > 0):
            macierz[i-1, i] = P1(i*h)
            macierz[i, i-1] = P3(i*h)
    macierz[n-1, n-1] = 1/h + 1/2
    "print(macierz)"

    "czy parzysty podzial ? "
    if (n%2 == 0):
        sr = (int)(n / 2) - 1
        macierz[sr,sr] = 2/h - 1/2
    else:
        l = (int)(n / 2) -1
        p = l + 1
        macierz[l,l] = 2/h - 1/8
        macierz[l,p] = -1/h - 7/8
        macierz[p,p] = 2/h - 1/8
        macierz[p,l] = -1/h + 5/8
    print("\n",  macierz)

    "metoda Simpsona"
    wektor_pr = np.zeros((n, 1))
    for i in range(1,n+1):
        wektor_pr[i-1] = f(2*i/n) * 4*h/3

    "zmiana na brzegach wektora prwaej strony"
    wektor_pr[0] = wektor_pr[0] - u0*(-1/h - 1/2)
    "wektor_pr[n-1] = wektor_pr[n-1] + g"
    wektor_pr[n-1] = (f(2) + 2 * f( 2- (h / 2))) * h / 6 + g

    print("\n",  wektor_pr)

    "rozwiazanie ukladu"
    rozwiazanie = la.solve(macierz, wektor_pr)

    wynik = np.zeros((n+1, 1))
    for i in range(n):
        wynik[i+1] = rozwiazanie[i]
    wynik[0] = u0
    "rysowanie"
    punkty = np.linspace(0.0, 2.0, n + 1)
    plt.plot(punkty, wynik)
    plt.show()
