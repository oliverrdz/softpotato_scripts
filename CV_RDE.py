import matplotlib.pyplot as plt
from softpotato.geometry import PlanarGeometry
from softpotato.reactions import ButlerVolmer, ElectrochemicalStep, Mechanism, Species
from softpotato.techniques import (
    CyclicVoltammetry,
    RotatingDiskElectrode,
    plot_sweep,
    simulate,
)

E0 = 0.0  # [V]
D = 1.0e-9  # [m^2/s]
c_bulk = 1.0  # [mol/m^3]
area = 1.96e-5  # [m^2]
nu = 1.0e-6  # [m^2/s]

ox = Species("Ox", diffusivity=D, bulk_concentration=c_bulk, standard_potential=E0, n_electrons=1)
red = Species("Red", diffusivity=D, bulk_concentration=0.0, standard_potential=E0, n_electrons=-1)
kinetics = ButlerVolmer(k0=5.0e-5, alpha=0.5, E0=E0, n=1)

mech = Mechanism(
    species=[ox, red],
    electrochemical_steps=[ElectrochemicalStep(ox, red, kinetics)],
)
geometry = PlanarGeometry(area=area)
cv_waveform = CyclicVoltammetry(E_start=0.35, E_vertex1=-0.45, E_vertex2=0.35, scan_rate=0.05)

rpms = [400.0, 900.0, 1600.0, 2500.0]
sweep_results = {}
for rpm in rpms:
    rde = RotatingDiskElectrode(
        waveform=cv_waveform,
        rpm=rpm,
        kinematic_viscosity=nu,
        model="nernst_layer",
    )
    sweep_results[f"{int(rpm)} rpm"] = simulate(rde, mech, geometry=geometry)

plot_sweep(sweep_results, legend_title="Rotation Speed")
plt.show()
