# Physical, numerical, etc constants

# Astronomical constants
from __future__ import annotations

# Earth constants
M_earth         = 5.972E24              # kg
R_earth         = 6.335439e6
R_core_earth    = 3485000.0             # m
M_core_earth    = 1.94E24               # kg
ocean_moles     = 7.68894973907177e+22 # moles of H2 (or H2O) in one present-day Earth ocean

# Physical constants
const_G = 6.67428e-11        #Gravitational constant (2006 measurements)
mol             = 6.02214076e+23        # mol definition
R_gas           = 8.31446261815324      # J K−1 mol−1

# Molar masses [kg mol-1]
molar_mass  = {
            "H"   : 0.001008,
            "C"   : 0.012011,
            "O"   : 0.015999,
            "N"   : 0.014007,
            "S"   : 0.03206,
            "He"  : 0.0040026,
            "H2O" : 0.01801528,
            "CO2" : 0.04401,
            "H2"  : 0.00201588,
            "CH4" : 0.01604,
            "CO"  : 0.02801,
            "N2"  : 0.028014,
            "O2"  : 0.031999,
            "SO2" : 0.064066,
            "H2S" : 0.0341,
            "S2"  : 0.0641,
            "NH3" : 0.017031,
        }

# Supported volatiles and elements
volatile_species = [ "H2O", "CO2", "H2", "CH4", "CO", "N2", "S2", "SO2", "H2S", "NH3"]
element_list     = [ "H", "O", "C", "N", "S" ]

# Plotting colours
dict_colors  = {
    "H2O": "#027FB1",
    "CO2": "#D24901",
    "H2" : "#008C01",
    "CH4": "#C720DD",
    "CO" : "#D1AC02",
    "N2" : "#870036",
    "S2" : "#FF8FA1",
    "SO2": "#00008B",
    "NH3": "#675200",
    "H2S": "#aaff22"
}
