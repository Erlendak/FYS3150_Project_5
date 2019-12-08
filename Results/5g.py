import matplotlib.pyplot as plt
import numpy as np
import sys

simulation_no_enrichments = np.loadtxt("simulation_no_enrichment.dat")
print("Reading simulation_no_enrichment.dat")

"""
Simulation with no enrichment.
Delta x ;  scaled by 120 Kilo meters

"""

#plt.plot(simulation_no_enrichments[8], label = r'$t_1$'+"Temperatur før likevekt.")
plt.plot(simulation_no_enrichments, label = "Temperatur ved likevekt." )
plt.legend()
plt.title("Simulasjon uten radiaktiv berikning ",size=17)
plt.ylabel("Temperatur ; Celsius",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
plt.grid()
plt.show()

"""
plt.imshow(simulation_no_enrichments,aspect='auto',cmap='hot_r')
plt.colorbar(label="Temperatur")
plt.title("Simulasjon uten radioaktiv berikning"+r'$\Delta X$'+" = ***",size=17)
plt.ylabel("Tid ; ",size=15)
plt.xlabel("Posisjon ; ",size=15)
plt.grid()
plt.show()
"""
