import matplotlib.pyplot as plt
import numpy as np
import sys
"""
result = []
forward_euler_dt01.dat = []
forward_euler_dt001.dat = []
backward_euler_dt01.dat = []
backward_euler_dt001.dat = []
crank_nicolson_dt01.dat = []
crank_nicolson_dt001.dat = []"""

"""
Mer eksempel kode.
A,B  = np.meshgrid(np.arange(9), np.zeros(12))
print(A.shape)
A = A.reshape(12, 3, 3)
print(A[0])
"""

"""
# result = np.array(result)
print(result.shape)
result  = result.astype(float).transpose()
xsize = 1000
ysize = 3
Nt = 3
result = result.reshape(xsize, ysize, Nt)
print(result.shape)"""


forward_euler_dt01 = np.loadtxt("forward_euler_dt01.dat")
forward_euler_dt001 = np.loadtxt("forward_euler_dt001.dat")
backward_euler_dt01 = np.loadtxt("backward_euler_dt01.dat")
backward_euler_dt001 = np.loadtxt("backward_euler_dt001.dat")
crank_nicolson_dt01 = np.loadtxt("crank_nicolson_dt01.dat")
crank_nicolson_dt001 = np.loadtxt("crank_nicolson_dt001.dat")


"""
Forward Euler.
Delta x ;  0.1

"""
plt.plot(forward_euler_dt01[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(forward_euler_dt01[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(forward_euler_dt01,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()


"""
Forward Euler.
Delta x ;  0.01

"""
plt.plot(forward_euler_dt001[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(forward_euler_dt001[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(forward_euler_dt001,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

"""
Backward Euler.
Delta x ;  0.1

"""
plt.plot(backward_euler_dt01[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(backward_euler_dt01[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(backward_euler_dt01,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()


"""
Backward Euler.
Delta x ;  0.01

"""
plt.plot(backward_euler_dt001[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(backward_euler_dt001[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(backward_euler_dt001,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()


"""
Crank Nicolson.
Delta x ;  0.1

"""
plt.plot(crank_nicolson_dt01[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(crank_nicolson_dt01[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Crank Nicolson "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(crank_nicolson_dt01,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Crank Nicolson "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()


"""
Crank Nicolson.
Delta x ;  0.01

"""
plt.plot(crank_nicolson_dt001[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(crank_nicolson_dt001[-1], label = r'$t_2$'+"Temperatur ved likevekt." )
plt.legend()
plt.title("Crank Nicolson "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Temperatur ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()

plt.imshow(crank_nicolson_dt001,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Crank Nicolson  "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()
