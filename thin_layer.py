import matplotlib.pyplot as plt
from softpotato.geometry import PlanarGeometry
from softpotato.grid import Uniform1DMesh
from softpotato.reactions import ButlerVolmer, Species
from softpotato.techniques import simulate, CyclicVoltammetry

ox = Species("O", diffusivity=1e-9, bulk_concentration=0)
red = Species("R", diffusivity=1e-9, bulk_concentration=1)
species = [ox, red]
geometry = PlanarGeometry(area=1e-5)
tech = CyclicVoltammetry(E_start=-0.4, E_vertex1=0.4, 
                         E_vertex2=-0.4, scan_rate=0.1)
kinetics = ButlerVolmer(k0=1e-3, alpha=0.5, E0=0, n=1)

grid = Uniform1DMesh(x_min=0, x_max=1e-5, n_points=101)
sim = simulate(tech, species, geometry, grid=grid, 
               kinetics=kinetics, bulk_bc="neumann")

plt.figure()
plt.plot(sim.E, sim.i*1e6)
plt.xlabel("$E$ / V")
plt.ylabel(r"$i$ / $\mu$A")
plt.grid(alpha=0.5)
plt.show()

k0_vals = [1e-4, 1e-5, 1e-6, 1e-7]
plt.figure()
for k0 in k0_vals:
    kinetics = ButlerVolmer(k0=k0, alpha=0.5, E0=0, n=1)
    sim = simulate(tech, species, geometry, grid=grid, kinetics=kinetics, bulk_bc="neumann")
    plt.plot(sim.E, sim.i*1e6, label=r"$k^0$ = " + f"{k0:.0e}" + r" cm s$^{-1}$")
plt.xlabel("$E$ / V")
plt.ylabel(r"$i$ / $\mu$A")
plt.grid(alpha=0.5)
plt.legend()
plt.show()



