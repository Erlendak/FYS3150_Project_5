import matplotlib.pyplot as plt
import numpy as np
import sys

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

plt.imshow(crank_nicolson_dt001,aspect='auto',cmap='hot_r', vmin = 0, vmax = 1 )
plt.colorbar(label="Temperatur")
plt.title("Crank Nicolson  "+r'$\Delta X$'+" = 0.01",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()


two_dimension_evolution_of_time = np.loadtxt("two_dimension_evolution_of_time.dat")
print("Reading two_dimension_evolution_of_time.dat")

two_dimension_evolution_of_time = two_dimension_evolution_of_time.reshape(4421, 10, 10)
"""print(two_dimension_evolution_of_time[0])
print("\n")
print(two_dimension_evolution_of_time[1])
""""""
print(two_dimension_evolution_of_time[2])
print(two_dimension_evolution_of_time[3])
print(two_dimension_evolution_of_time[-1])
print(two_dimension_evolution_of_time[-2])
"""

def show_temperature_at(_t):
    plt.imshow(two_dimension_evolution_of_time[_t],aspect='auto',cmap='hot_r',extent=[0,10,0,10])
    plt.colorbar(label="Temperatur")
    plt.title("Temperatur i et kvadrat, ved tiden t = {}.\n".format(_t*0.000123457)+r'$\Delta X, \Delta Y$'+" = 0.1111...    "+r'$\Delta t$'+" = 0.000123457...",size=17)
    plt.ylabel("Posisjon Y-akse ; ",size=15)
    plt.xlabel("Posisjon X-akse ; ",size=15)
    #plt.xlim([0,9])
    #plt.ylim([0,9])
    plt.yticks(np.arange(0, 10, 1))
    plt.xticks(np.arange(0, 10, 1))
    plt.grid()#which = 'minor',linewidth=1)

    plt.show()

show_temperature_at(0)
show_temperature_at(1)
show_temperature_at(2)
show_temperature_at(3)
show_temperature_at(49)
show_temperature_at(499)
show_temperature_at(4419)
show_temperature_at(4420)

two_dimension = np.loadtxt("two_dimension.dat")
print("Reading two_dimension.dat")
two_dimension = two_dimension.reshape(2, 10, 10)

plt.imshow(two_dimension[1],aspect='auto',cmap='hot_r',extent=[0,10,0,10])
plt.colorbar(label="Temperatur")
plt.title("Temperatur i et kvadrat, ved tiden t = {}.\n".format("eqvilibrium")+r'$\Delta X, \Delta Y$'+" = 0.1111...    "+r'$\Delta t$'+" = 0.000123457...",size=17)
plt.ylabel("Posisjon Y-akse ; ",size=15)
plt.xlabel("Posisjon X-akse ; ",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(0, 10, 1))
plt.xticks(np.arange(0, 10, 1))
plt.grid()#which = 'minor',linewidth=1)
plt.show()
