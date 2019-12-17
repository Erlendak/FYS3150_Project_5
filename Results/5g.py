import matplotlib.pyplot as plt
import numpy as np
import sys

simulation_no_enrichments = np.loadtxt("simulation_no_enrichment.dat")
print("Reading simulation_no_enrichment.dat")

"""
Simulation with no enrichment.
Delta x ;  1 Kilometer

"""

plt.plot(np.linspace(0,120,200),np.linspace(8,1300,200), label = "Temperatur uten varmeproduksjon")
plt.plot(simulation_no_enrichments, label = "Temperatur med varmeproduksjon ved likevekt." )
plt.legend()
plt.title("Simulasjon med varmeproduksjon.\nUten radioaktiv berikning.",size=17)
plt.ylabel("Temperatur ; Celsius",size=15)
plt.xlabel("Dybde ; Kilometer",size=15)
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

print("Reading simulation_with_enrichment.dat")
simulation_with_enrichments = np.loadtxt("simulation_with_enrichment.dat")
simulation_with_enrichments = simulation_with_enrichments.reshape(2, 120, 150)
"""
Simulation with enrichment.
Lenght X ; 120 Kilometer
Lenght Y ; 150 Kilometer
Delta x ;  1 Kilometer
Delta y ;  1 Kilometer

"""
plt.subplot(1,2,1)
plt.imshow(simulation_with_enrichments[0],aspect='auto',cmap='hot_r',extent=[0,150,120,0])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Før radioaktiv berikning med kun varmeproduksjon\nVed oppnådd equilibrium stadie.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(0, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)

plt.subplot(1,2,2)
plt.imshow(simulation_with_enrichments[1],aspect='auto',cmap='hot_r',extent=[0,150,120,0])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Varmeproduksjon med konstant radioaktiv berikning\nVed oppnådd equilibrium stadie.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(0, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)
plt.show()

plt.subplot(1,2,1)
plt.imshow(simulation_with_enrichments[0][0:150][60:120],aspect='auto',cmap='hot_r',extent=[0,150,120,60])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Mantelen før radiaktiv berikning\nInitial betingelsene til simasjonen med radioaktiv berikning",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(60, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)

plt.subplot(1,2,2)
plt.imshow(simulation_with_enrichments[1][0:150][60:120],aspect='auto',cmap='hot_r',extent=[0,150,120,60])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Mantelen med radioaktiv berikning\nVed ekvilibrium.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(60, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)
plt.show()
#print(simulation_with_enrichments[1])

print("Reading simulation_implemented_enrichment_decay.dat")
simulation_implemented_enrichment_decay = np.loadtxt("simulation_implemented_enrichment_decay.dat")
simulation_implemented_enrichment_decay = simulation_implemented_enrichment_decay.reshape(10, 120, 150)

plt.subplot(1,2,1)
plt.imshow(simulation_with_enrichments[1],aspect='auto',cmap='hot_r',extent=[0,150,120,0])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Temperatur distrubisjonen i lithospheren\nFor et Giga år siden.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(0, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)

plt.subplot(1,2,2)
plt.imshow(simulation_implemented_enrichment_decay[9],aspect='auto',cmap='hot_r',extent=[0,150,120,0])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Temperatur distrubisjonen i lithospheren\nSimulert fra et Giga år siden til i dag.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(0, 120, 20))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)
plt.show()

plt.subplot(1,2,1)
plt.imshow(simulation_with_enrichments[1][0:150][60:120],aspect='auto',cmap='hot_r',extent=[0,150,120,60])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Temperatur distrubisjonen i mantelen\nFor et Giga år siden.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(60, 120, 10))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)

plt.subplot(1,2,2)
plt.imshow(simulation_implemented_enrichment_decay[9][0:150][60:120],aspect='auto',cmap='hot_r',extent=[0,150,120,60])
plt.colorbar(label="Temperatur; Celsius")
plt.title("Temperatur distrubisjonen i mantelen\nSimulert fra et Giga år siden til i dag.",size=17)
plt.ylabel("Dybde ; Kilometer",size=15)
plt.xlabel("Posisjon ; Kilometer",size=15)
#plt.xlim([0,9])
#plt.ylim([0,9])
plt.yticks(np.arange(60, 120, 10))
plt.xticks(np.arange(0, 150, 20))
plt.grid()#which = 'minor',linewidth=1)
plt.show()
