import matplotlib.pyplot as plt
from softpotato.geometry import PlanarGeometry
from softpotato.reactions import (
    ButlerVolmer,
    ElectrochemicalStep,
    FirstOrderReaction,
    Mechanism,
    NernstianEquilibrium,
    Species,
)
from softpotato.techniques import SampledCurrentVoltammetry, plot_scv, simulate

# Physical and transport parameters (SI units)
D = 1e-9          # Diffusion coefficient [m^2/s] (1e-5 cm^2/s)
c_bulk = 1.0      # Bulk concentration of Ox [mol/m^3] (1 mM)
area = 1e-4       # Planar electrode area [m^2] (1 cm^2)
E0 = 0.0          # Formal reduction potential [V]
n = 1             # Number of electrons transferred
k_chem = 2.0      # Homogeneous follow-up rate constant [s^-1]

# Define chemical species
ox = Species("Ox", diffusivity=D, bulk_concentration=c_bulk, standard_potential=E0, n_electrons=n)
red = Species("Red", diffusivity=D, bulk_concentration=0.0, standard_potential=E0, n_electrons=-n)
prod = Species("P", diffusivity=D, bulk_concentration=0.0)

# Define the EC mechanism:
# Use NernstianEquilibrium for reversible E-step (ErCi), or swap in ButlerVolmer for finite kinetics
e_kinetics = NernstianEquilibrium(E0=E0, n=n)
# e_kinetics = ButlerVolmer(k0=1e-4, alpha=0.5, E0=E0, n=n)

ec_mechanism = Mechanism(
    species=[ox, red, prod],
    reactions=[FirstOrderReaction(red, prod, k=k_chem)],
    electrochemical_steps=[ElectrochemicalStep(ox, red, e_kinetics)],
    name="EC_Mechanism",
)

# Configure Sampled-Current Voltammetry (potential steps & sampling times tau)
scv = SampledCurrentVoltammetry(
    E_start=0.25,
    E_end=-0.25,
    dE=0.025,
    duration=2.0,
    sample_times=[0.1, 0.25, 0.5, 1.0, 2.0],
)

# Run event-aware PDE simulation
result = simulate(scv, ec_mechanism, geometry=PlanarGeometry(area=area))

# Extract sampled voltammogram arrays: E_steps (N_steps,) and i_sampled (N_tau, N_steps)
E_steps, i_sampled = result.sample_currents()

# Plot sampled-current voltammograms across all sampling times tau
fig, ax = plot_scv(result)
ax.set_title(f"Sampled-Current Voltammogram — EC Mechanism ($k_c = {k_chem}\\,\\mathrm{{s}}^{{-1}}$)")
plt.tight_layout()
plt.show()
