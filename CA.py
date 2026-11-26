from softpotato.geometry import PlanarGeometry
from softpotato.reactions import NernstianEquilibrium, Species
from softpotato.techniques import Chronoamperometry, simulate

ox = Species("O", diffusivity=1e-9, bulk_concentration=1)
red = Species("R", diffusivity=1e-9, bulk_concentration=0.0)
species = [ox, red]
kinetics = NernstianEquilibrium(E0=0, n=1)
geometry = PlanarGeometry(area=1e-4)
technique = Chronoamperometry(E_step=-0.5, duration=1)

sim = simulate(technique, species, geometry, kinetics=kinetics)
sim.plot_summary()
