# Physical, numerical, etc constants

# Astronomical constants
from __future__ import annotations

# Earth constants
M_earth = 5.972e24  # kg
R_earth = 6.335439e6
R_core_earth = 3485000.0  # m
M_core_earth = 1.94e24  # kg
# Moles of H2 (or H2O) in one present-day Earth ocean.
# Source: Clark, W. C. (1982), Carbon Dioxide Review, p. 469, Oxford Univ. Press, New York.
ocean_moles = 7.68894973907177e22

# Physical constants
const_G = 6.67428e-11  # Gravitational constant (2006 measurements)
mol = 6.02214076e23  # mol definition
R_gas = 8.31446261815324  # J K−1 mol−1

# Molar masses [kg mol-1]
molar_mass = {
    'H': 0.001008,
    'C': 0.012011,
    'O': 0.015999,
    'N': 0.014007,
    'S': 0.03206,
    'He': 0.0040026,
    'Ne': 0.0201797,
    'Ar': 0.039948,
    'Kr': 0.083798,
    'Xe': 0.131293,
    'H2O': 0.01801528,
    'CO2': 0.04401,
    'H2': 0.00201588,
    'CH4': 0.01604,
    'CO': 0.02801,
    'N2': 0.028014,
    'O2': 0.031999,
    'SO2': 0.064066,
    'H2S': 0.0341,
    'S2': 0.0641,
    'NH3': 0.017031,
}

# Noble gases. Each is monatomic and chemically inert: its gas species and
# its element are the same entity, it takes no part in the CHNOS reaction
# network, and it partitions between melt and atmosphere by Henry's law
# alone. Ordered by atomic number so the active-species vector the solver
# builds is deterministic. They are part of `element_list` (the full set of
# supported elements) but not of `volatile_species` (the reaction network);
# a noble gas with no budget is inactive and contributes nothing.
noble_gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe']

# Supported reaction-network volatiles.
volatile_species = ['H2O', 'CO2', 'O2', 'H2', 'CH4', 'CO', 'N2', 'S2', 'SO2', 'H2S', 'NH3']

# The reacting elements of the CHNOS network. The solver's element-residual and
# element-mass loops iterate this subset; the inert noble gases partition by
# Henry's law and are handled on their own active-gas path.
element_list_chnos = ['H', 'O', 'C', 'N', 'S']

# All supported elements: the reacting CHNOS elements plus the inert noble
# gases. This is the complete element registry.
element_list = element_list_chnos + noble_gases

# Plotting colours
dict_colors = {
    'H2O': '#027FB1',
    'CO2': '#D24901',
    'O2': '#00dd00',
    'H2': '#008C01',
    'CH4': '#C720DD',
    'CO': '#D1AC02',
    'N2': '#870036',
    'S2': '#FF8FA1',
    'SO2': '#00008B',
    'NH3': '#675200',
    'H2S': '#aaff22',
    'He': '#7f7f7f',
    'Ne': '#e6550d',
    'Ar': '#3182bd',
    'Kr': '#31a354',
    'Xe': '#756bb1',
}
