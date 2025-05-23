import matplotlib.pyplot as plt

g = 32.17
Cd = 0.75
A = 78.54/144
a = 0
v = 0
Alt = 1
M_tot = (155.69 + 10 + 12.65*2 + 31.64*2)/g
Mdot = 4.43/g
rho = 4.82331 * 10**(-13) * Alt**2 - 6.37675 * 10**(-8) * Alt + 0.002358692
F_t = 1000
F_g = M_tot*g
F_d = 0.5 * rho * v**2 * Cd * A
F_net = F_t - F_g - F_d
FlightTime = 0

alist = [1]
vlist = [0]
Altlist = [0]
Mlist = [M_tot*g]
Flist = [F_net]
rholist = [0.002358692]
TimeList = [0]

while FlightTime < 50:
    print(f"F_net: {F_net} | M_tot: {M_tot} | FlightTime: {FlightTime} | Vel: {v} | Alt: {Alt}")
    a = F_net/M_tot
    v = v + a
    Alt = Alt + v



    if FlightTime > 19:
        F_t = 0
    else:
        M_tot = M_tot - Mdot
        F_g = M_tot*g

    rho = 4.82331 * 10**(-13) * Alt**2 - 6.37675 * 10**(-8) * Alt + 0.002358692

    if v >= 0:
        F_d = 0.5 * rho * v**2 * Cd * A
    else: 
        F_d = -0.5 * rho * v**2 * Cd * A

    F_net = F_t - F_g - F_d

    FlightTime += 1

    alist.append(a)
    vlist.append(v)
    Altlist.append(Alt)
    Mlist.append(M_tot*g)
    Flist.append(F_net)
    rholist.append(rho)
    TimeList.append(FlightTime)

plt.figure()
plt.plot(TimeList, alist, label="a")
plt.xlabel("time")
plt.ylabel("a")
plt.legend()
plt.show(block=False)
plt.figure()
plt.plot(TimeList, vlist, label="v")
plt.xlabel("time")
plt.ylabel("v")
plt.legend()
plt.show(block=False)
plt.figure()
plt.plot(TimeList, Altlist, label="Alt")
plt.xlabel("time")
plt.ylabel("Alt")
plt.legend()
plt.show(block=False)
plt.figure()
plt.plot(TimeList, Mlist, label="m")
plt.xlabel("time")
plt.ylabel("m")
plt.legend()
plt.show(block=False)
plt.figure()
plt.plot(TimeList, Flist, label="F")
plt.xlabel("time")
plt.ylabel("F")
plt.legend()
plt.show(block=False)
plt.figure()
plt.plot(TimeList, rholist, label="rho")
plt.xlabel("time")
plt.ylabel("rho")
plt.legend()
plt.show(block=False)
plt.show()


