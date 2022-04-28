import numpy as np
import matplotlib.pyplot as plt

real_avg_mag = []
avg_mag_sqrd = []

maxT = 2
l = 100
T = np.linspace(0, (l-1)*maxT/l, l)

with open("maggy.txt", "r") as infile:
    for i in range(l):
        line = infile.readline().split(",")
        real_avg_mag.append(float(line[0]))
        avg_mag_sqrd.append(float(line[1]))

ram = np.array(real_avg_mag)
ams = np.array(avg_mag_sqrd)

plt.figure()
plt.plot(T, ram)
plt.title("Real part of average magnetization")

plt.figure()
plt.plot(T, ams)
plt.title("Absolute magnetization squared")

plt.show()
