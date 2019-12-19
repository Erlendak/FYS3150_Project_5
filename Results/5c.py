import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import numpy as np
import sys


dt = 0.00005
t1 =int( 0.1/dt)


forward_euler_dt01 = np.loadtxt("forward_euler_dt01.dat")
print("Reading forward_euler_dt01.dat")
forward_euler_dt001 = np.loadtxt("forward_euler_dt001.dat")
print("Reading forward_euler_dt001.dat")
backward_euler_dt01 = np.loadtxt("backward_euler_dt01.dat")
print("Reading backward_euler_dt01.dat")
backward_euler_dt001 = np.loadtxt("backward_euler_dt001.dat")
print("Reading backward_euler_dt001.dat")
crank_nicolson_dt01 = np.loadtxt("crank_nicolson_dt01.dat")
print("Reading crank_nicolson_dt01.dat")
crank_nicolson_dt001 = np.loadtxt("crank_nicolson_dt001.dat")
print("Reading crank_nicolson_dt001.dat")


"""
Analytical solution

"""

for nx in [11,101]:
    x = np.linspace(0,1,nx)

    res = np.zeros(nx)
    t=0.10
    for k in range(1,nx):
        res += 2*(-1)**k/(k*np.pi)*np.exp(-(k*np.pi)**2*t)*np.sin(np.pi*k*x)
    res += x

for nx2 in [11,101]:
    x2 = np.linspace(0,1,nx2)

    res2 = np.zeros(nx2)
    t2= 1
    for k in range(1,nx2):
        res2 += 2*(-1)**k/(k*np.pi)*np.exp(-(k*np.pi)**2*t2)*np.sin(np.pi*k*x2)
    res2 += x2

plt.subplot(2,2,1)
plt.plot(x,res,label= r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(x,res2,label=r'$t_2$'+" = 1.0 s,  Ved likevekt.")
plt.title("Analytisk løsning",size=17)
plt.ylabel("u ( x, t ) ;",size=15)
plt.xlabel("x ; ",size=15)
plt.legend()
plt.grid()


"""
Forward Euler.
Delta x ;  0.1

"""

plt.subplot(2,2,2)
plt.plot(forward_euler_dt01[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(forward_euler_dt01[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()


"""
Backward Euler.
Delta x ;  0.1

"""

plt.subplot(2,2,3)
plt.plot(backward_euler_dt01[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(backward_euler_dt01[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()


"""
Crank Nicolson.
Delta x ;  0.1

"""

plt.subplot(2,2,4)
plt.plot(crank_nicolson_dt01[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(crank_nicolson_dt01[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Crank Nicolson "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()
plt.subplots_adjust(hspace=0.3)
plt.show()




"""
Analytical solution.

"""

plt.subplot(2,2,1)
plt.plot(x,res,label= r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(x,res2,label=r'$t_2$'+" = 1.0 s,  Ved likevekt.")
plt.title("Analytisk løsning",size=17)
plt.ylabel("u ( x, t ) ;",size=15)
plt.xlabel("x ; ",size=15)
plt.legend()
plt.grid()


"""
Forward Euler.
Delta x ;  0.01

"""

plt.subplot(2,2,2)
plt.plot(forward_euler_dt001[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(forward_euler_dt001[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Forward Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()


"""
Backward Euler.
Delta x ;  0.01

"""

plt.subplot(2,2,3)
plt.plot(backward_euler_dt001[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(backward_euler_dt001[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Backwards Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()


"""
Crank Nicolson.
Delta x ;  0.01

"""

plt.subplot(2,2,4)
plt.plot(crank_nicolson_dt001[t1], label = r'$t_1$'+" = 0.1 s,  Før likevekt.")
plt.plot(crank_nicolson_dt001[-1], label = r'$t_2$'+" = 1.0 s,  Ved likevekt." )
plt.legend()
plt.title("Crank Nicolson "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("u ( x, t ) ; ",size=15)
plt.xlabel("x ; ",size=15)
plt.grid()
plt.subplots_adjust(hspace=0.3)
plt.show()




"""
Imshow

"""


"""
Forward Euler.
Delta x ;  0.1

"""

plt.subplot(3,2,1)
plt.imshow(forward_euler_dt01,aspect='auto',cmap='hot_r',extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nForward Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()


"""
Forward Euler.
Delta x ;  0.01

"""

plt.subplot(3,2,2)
plt.imshow(forward_euler_dt001,aspect='auto',cmap='hot_r',extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nForward Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()


"""
Backward Euler.
Delta x ;  0.1

"""

plt.subplot(3,2,3)
plt.imshow(backward_euler_dt01,aspect='auto',cmap='hot_r',extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nBackwards Euler "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()


"""
Backward Euler.
Delta x ;  0.01

"""

plt.subplot(3,2,4)
plt.imshow(backward_euler_dt001,aspect='auto',cmap='hot_r',extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nBackwards Euler "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()


"""
Crank Nicolson.
Delta x ;  0.1

"""

plt.subplot(3,2,5)
plt.imshow(crank_nicolson_dt01,aspect='auto',cmap='hot_r',extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nCrank Nicolson "+r'$\Delta X$'+" = 0.1",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()


"""
Crank Nicolson.
Delta x ;  0.01

"""

plt.subplot(3,2,6)
plt.imshow(crank_nicolson_dt001,aspect='auto',cmap='hot_r', vmin = 0, vmax = 1,extent=[0,10,1,0])
plt.colorbar(label=" u ( x , t ) ;")
plt.title("u ( x , t ) \nCrank Nicolson  "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel(" t ; ",size=15)
plt.xlabel(" x ; ",size=15)
plt.grid()
plt.subplots_adjust(hspace=0.3)
plt.show()



"""
2. Dimentional example.

"""

two_dimension_evolution_of_time = np.loadtxt("two_dimension_evolution_of_time.dat")
print("Reading two_dimension_evolution_of_time.dat")

two_dimension_evolution_of_time = two_dimension_evolution_of_time.reshape(4421, 10, 10)

n = 10
dx = 1.0/(n-1)   #// Set delta x for the simulation.
dt = 0.01*dx*dx  #// set appropriate amount of time step.

def show_temperature_at(_t):
    plt.imshow(rotate(two_dimension_evolution_of_time[_t],-90),aspect='equal',cmap='hot_r',extent=[0,10,0,10])
    plt.colorbar(label="U ( X, Y )")
    plt.title("U ( X, Y, t )\nVed tiden t = {0:.2f} s.\n".format(_t*dt),size=17)
    plt.ylabel("Posisjon Y-akse ; Y ",size=15)
    plt.xlabel("Posisjon X-akse ; X ",size=15)
    plt.yticks(np.arange(0, 10, 1))
    plt.xticks(np.arange(0, 10, 1))
    plt.grid()#which = 'minor',linewidth=1)


plt.subplot(2,2,1)
show_temperature_at(0)

plt.subplot(2,2,2)
show_temperature_at(49)

plt.subplot(2,2,3)
show_temperature_at(499)

plt.subplot(2,2,4)
show_temperature_at(4420)

plt.show()


"""
Check equilibrium state should be same as last step in evolution of time
""" """

two_dimension = np.loadtxt("two_dimension.dat")
print("Reading two_dimension.dat")
two_dimension = two_dimension.reshape(2, 10, 10)

plt.imshow(two_dimension[1],aspect='auto',cmap='hot_r',extent=[0,10,0,10])
plt.colorbar(label="U ( X, Y )")
plt.title("U ( X, Y, t )\nVed equilibrium.\n",size=17)
plt.ylabel("Posisjon Y-akse ; Y ",size=15)
plt.xlabel("Posisjon X-akse ; X ",size=15)
plt.yticks(np.arange(0, 10, 1))
plt.xticks(np.arange(0, 10, 1))
plt.grid()
plt.show()
"""
