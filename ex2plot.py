import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes as zia

maxT = 2
l = 100
T = np.linspace(0, (l-1)*maxT/l, l)

def get_mags(filename):
    with open(filename, "r") as infile:
        real_avg_mag = []
        avg_mag_sqrd = []
        avg_mag_fth = []

        for i in range(l):
            line = infile.readline().split(",")
            real_avg_mag.append(float(line[0]))
            avg_mag_sqrd.append(float(line[1]))
            avg_mag_fth.append(float(line[2]))

    ram = np.array(real_avg_mag)
    ams = np.array(avg_mag_sqrd)
    amf = np.array(avg_mag_fth)
    gamma = amf/ams**2

    return ram, ams, gamma

ram8, ams8, gamma8 = get_mags("maggy8.txt")
ram16, ams16, gamma16 = get_mags("maggy16.txt")
ram32, ams32, gamma32 = get_mags("maggy32.txt")

r_ram8, r_ams8, r_gamma8 = get_mags("r_maggy8.txt")
r_ram16, r_ams16, r_gamma16 = get_mags("r_maggy16.txt")
r_ram32, r_ams32, r_gamma32 = get_mags("r_maggy32.txt")

plt.figure()
plt.plot(T, ram16)
plt.xlabel("$T/J$")
plt.ylabel("$\Re{\langle m \\rangle}$")
plt.title("Real part of average magnetization\n$L = 16$, Case A")
plt.savefig("ram.png")

plt.figure()
plt.plot(T, r_ram16)
plt.xlabel("$T/J$")
plt.ylabel("$\Re{\langle m \\rangle}$")
plt.title("Real part of average magnetization\n$L = 16$, Case B")
plt.savefig("r_ram.png")

plt.figure()
plt.plot(T, ams16)
plt.xlabel("$T/J$")
plt.ylabel("$\langle |m|^2 \\rangle$")
plt.title("Absolute magnetization squared\n$L = 16$")
plt.savefig("ams.png")

fig = plt.figure()
#ax = plt.gca()
plt.plot(T, gamma8, label = "$L = 8$")
plt.plot(T, gamma16, label = "$L = 16$")
plt.plot(T, gamma32, label = "$L = 32$")
plt.axvline(1/(np.log(1+np.sqrt(3))), c="r", ls="--", label = "$T_c/J$",\
        alpha = 0.4)
plt.legend()

#axins = zia(ax, zoom = 1, loc="lower right")
#axins.set_xlim(0.95, 1.05)
#axins.set_ylim(1.05, 1.2)

plt.xlabel("$T/J$")
plt.ylabel("$\Gamma(T/J)$")
plt.title("Gamma as a function of temperature")
plt.savefig("gamma.png")

'''
fig = plt.figure()
plt.plot(T, gamma8, label = "$L = 8$")
plt.plot(T, gamma16, label = "$L = 16$")
plt.plot(T, gamma32, label = "$L = 32$")
plt.axvline(1/(np.log(1+np.sqrt(3))), c="r", ls="--", label = "$T_c/J$",\
        alpha = 0.4)
plt.legend()
plt.xlabel("$T/J$")
plt.ylabel("$\Gamma(T/J)$")
plt.title("Gamma as a function of temperature")
#plt.savefig("r_gamma.png")
'''

plt.show()
