from softpotato.geometry import PlanarGeometry
from softpotato.reactions import NernstianEquilibrium, Species
from softpotato.techniques import CyclicVoltammetry, simulate

import matplotlib.pyplot as plt

ox = Species("O", diffusivity=1e-9, bulk_concentration=0)
red = Species("R", diffusivity=1e-9, bulk_concentration=1)
species = [ox, red]
kinetics = NernstianEquilibrium(E0=0, n=1)
geometry = PlanarGeometry(area=1e-4)
technique = CyclicVoltammetry(E_start=-0.4, E_vertex1=0.4, E_vertex2=-0.4, 
                              scan_rate=0.1)

sim = simulate(technique, species, geometry, kinetics=kinetics)
sim.plot_summary()
plt.show()
