import numpy as np
import matplotlib.pyplot as plt

mag = []
amag2 = []

l = 100
maxtemp = 2

with open("1Dmaggy.txt", "r") as infile:
    for i in range(100):

        line = infile.readline().split(",")
        mag.append(float(line[0]))
        amag2.append(float(line[1]))

mag = np.array(mag)
amag2 = np.array(amag2)

t = np.linspace(0, (l-1)*maxtemp/l, l)

plt.figure()
plt.title("Real part of average magnetization")
plt.plot(t, mag)

plt.figure()
plt.title("Absolute square of the average magnetization")
plt.plot(t, amag2)

plt.show()
