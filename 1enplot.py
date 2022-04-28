import numpy as np
import matplotlib.pyplot as plt

def U(beta, J, L):
    e = np.exp(beta*J)
    return (-L*J*e * (2*(e-1)**(L-1) + (e+2)**(L-1)))/(2*(e-1)**L + (e+2)**L)

def free(beta, J, L):
    e = np.exp(beta*J)
    return 1/beta * np.log(2*(e-1)**L + (e+2)**L)

beta = np.linspace(0, 5, 1000)

for l in [4, 8, 16, 32, 64]:
    plt.plot(beta, U(beta, 1, l)/l, label = f"Energy, L = {l}")

plt.legend()
plt.show()

plt.plot(beta, free(beta, 1, 8))
plt.show()
