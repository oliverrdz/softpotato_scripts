from softpotato.geometry import PlanarGeometry
from softpotato.reactions import NernstianEquilibrium, Species
from softpotato.techniques import simulate

import matplotlib.pyplot as plt 

ox = Species("O", diffusivity=1e-9, bulk_concentration=0)
red = Species("R", diffusivity=1e-9, bulk_concentration=1)
species = [ox, red]
kinetics = NernstianEquilibrium(E0=0, n=1)
geometry = PlanarGeometry(area=1e-4)

# Chronoamperometry
from softpotato.techniques import Chronoamperometry
tech = Chronoamperometry(E_step=0.5, duration=1)
sim = simulate(tech, species, geometry, kinetics=kinetics)
sim.plot()

# Cyclic Voltammetry
from softpotato.techniques import CyclicVoltammetry
tech = CyclicVoltammetry(E_start=-0.4, E_vertex1=0.4, 
                         E_vertex2=-0.4, scan_rate=0.1)
sim = simulate(tech, species, geometry, kinetics=kinetics)
sim.plot()

# Sampled Current Voltammetry
from softpotato.techniques import SampledCurrentVoltammetry
tech = SampledCurrentVoltammetry(E_start=-0.4, E_end=0.4, dE=0.01, 
                                 duration=1, sample_times=[0.1, 0.5, 1])
sim = simulate(tech, species, geometry, kinetics=kinetics)
sim.plot()

# Rotating Disc Electrode
from softpotato.techniques import RotatingDiskElectrode
tech = RotatingDiskElectrode(E_start=-0.4, E_end=0.4, scan_rate=0.01, rpm=1000)
sim = simulate(tech, species, geometry, kinetics=kinetics)
sim.plot()

plt.show()

