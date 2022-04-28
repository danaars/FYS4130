import numpy as np
import matplotlib.pyplot as plt

def C(r, L = 16, beta = 6, J = 1):  # Analytical expression
    e = np.exp(beta*J)
    return ((e-1)**L + (e-1)**r * (e+2)**(L-r) + \
            (e-1)**(L-r) * (e+2)**r) / (2*(e-1)**L + (e+2)**L)

r = np.arange(0, 16)
#r = np.linspace(0, 16, 1000)

hgT = []
lwT = []

with open("corr1D_0.500000.txt", "r") as infile:
    for i in range(16):
        hgT.append(float(infile.readline()))

with open("corr1D_0.250000.txt", "r") as infile:
    for i in range(16):
        lwT.append(float(infile.readline()))

plt.plot(r, C(r, beta = 2), label = "Analytical")
plt.plot(r, hgT, label = "Numerical")
plt.title("High Temperature case")
plt.legend()
plt.show()

plt.plot(r, C(r, beta = 4), label = "Analytical")
plt.plot(r, lwT, label = "Numerical")
plt.title("Low Temperature case")
plt.legend()
plt.show()
