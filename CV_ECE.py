import matplotlib.pyplot as plt
from softpotato.geometry import PlanarGeometry
from softpotato.reactions import (
    ButlerVolmer,
    ElectrochemicalStep,
    FirstOrderReaction,
    Mechanism,
    Species,
)
from softpotato.techniques import CyclicVoltammetry, plot_result, simulate

# 1. Define chemical species in SI units
# D in m^2/s (1e-9 m^2/s = 1e-5 cm^2/s), bulk_concentration in mol/m^3 (1.0 mol/m^3 = 1 mM)
D = 1e-9
species_A = Species("A", diffusivity=D, bulk_concentration=1.0, n_electrons=1)
species_B = Species("B", diffusivity=D, bulk_concentration=0.0, n_electrons=-1)
species_C = Species("C", diffusivity=D, bulk_concentration=0.0, n_electrons=1)
species_D = Species("D", diffusivity=D, bulk_concentration=0.0, n_electrons=-1)

# 2. Define heterogeneous (E1, E2) and homogeneous (C) reaction steps
# E1: A + e- <-> B at E0_1 = 0.0 V
e1_kinetics = ButlerVolmer(k0=1e-3, alpha=0.5, E0=0.0, n=1)  # k0 in m/s

# C: B -> C with first-order rate constant k_c = 5.0 s^-1
chem_step = FirstOrderReaction(species_B, species_C, k=5.0)

# E2: C + e- <-> D at E0_2 = -0.25 V
e2_kinetics = ButlerVolmer(k0=1e-3, alpha=0.5, E0=-0.25, n=1)

ece_mechanism = Mechanism(
    species=[species_A, species_B, species_C, species_D],
    reactions=[chem_step],
    electrochemical_steps=[
        ElectrochemicalStep(species_A, species_B, e1_kinetics),
        ElectrochemicalStep(species_C, species_D, e2_kinetics),
    ],
    name="ECE_Mechanism",
)

# 3. Configure electrode geometry and cyclic voltammetry waveform
geometry = PlanarGeometry(area=1e-4)  # m^2 (1 cm^2)
cv = CyclicVoltammetry(
    E_start=0.3,
    E_vertex1=-0.5,
    E_vertex2=0.3,
    scan_rate=0.1,  # V/s
)

# 4. Run the event-aware PDE simulation and plot the voltammogram
result = simulate(cv, ece_mechanism, geometry=geometry)

fig, ax = plot_result(result)
ax.set_title("ECE Mechanism Cyclic Voltammogram")
plt.show()
