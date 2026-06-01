"""Authoritative elemental inventories used by the cross-backend
comparison harness.

Earth bulk-silicate-Earth (BSE) H / C / N / S inventory is taken from
Krijt et al. (2023), Protostars and Planets VII, Tables 1 and 2 ("BSE"
totals row).

Oxygen is a special case. Krijt et al. tabulate redox-active O
following the Evans (2006) convention (mass of O required to move the
silicate Earth to Fe(II)O reference state), which is dominated by the
mantle FeO / Fe2O3 redox imbalance, not by volatile O. The
authoritative-O entry point treats O as the volatile budget (the O
atoms residing in atmospheric H2O / CO2 / SO2 / O2 and dissolved as
the same species). The two definitions are inequivalent and not
interconvertible without a chemistry calculation.

We therefore derive the canonical Earth volatile-O reference from the
Krijt BSE H/C/N/S budget by running the buffered-mode solver at
T_magma = 2000 K and Delta-IW = +3.5 (the Sossi 2020 estimate of
Earth's modern upper-mantle fO2). The resulting O_kg_total is the
self-consistent volatile O for Earth at that thermodynamic state. The
value is stored as `EARTH_VOLATILE_O_REF_KG` and used as the 1x point
of the Fig 3 grid sweep.

The numbers here are bulk-silicate-Earth, not bulk-Earth. The core is
excluded.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Inventory:
    """Elemental inventory in kg, with citation provenance."""

    name: str
    H: float
    C: float
    N: float
    O: float  # noqa: E741  O is the oxygen-element field name across CALLIOPE
    S: float
    citation: str
    notes: str = ''

    def asdict(self) -> dict:
        """Return as `target_d` shape for `equilibrium_atmosphere_authoritative_O`."""
        return {'H': self.H, 'C': self.C, 'N': self.N, 'O': self.O, 'S': self.S}


# Krijt et al. (2023) PPVII Tables 1+2 BSE H/C/N/S in kg. Used as the
# authoritative H/C/N/S inventory for the Earth fiducial; oxygen is
# derived self-consistently from a chemistry call (see module docstring).
EARTH_HCNS_KRIJT23 = {
    'H': 5.6e20,
    'C': 3.1e21,
    'N': 3.7e19,
    'S': 1.0e21,
}


# Earth's volatile O at the Sossi 2020 Delta-IW = +3.5, T = 2000 K state.
# Computed by `derive_earth_volatile_O()` from CALLIOPE's buffered-mode
# solver with EARTH_HCNS_KRIJT23 as the H/C/N/S target and the current
# default Fischer 2011 IW buffer. Hard-coded here so the harness does
# not have to recompute it on every invocation; the provenance script
# `derive_earth_volatile_O()` re-derives it on demand to confirm the
# constant has not drifted. The legacy O'Neill 2002 buffer gives a
# slightly different value (~1.241e22 kg).
EARTH_VOLATILE_O_REF_KG = 1.260e22


EARTH_BSE_KRIJT23 = Inventory(
    name='Earth BSE (volatile O at IW+3.5)',
    H=EARTH_HCNS_KRIJT23['H'],
    C=EARTH_HCNS_KRIJT23['C'],
    N=EARTH_HCNS_KRIJT23['N'],
    O=EARTH_VOLATILE_O_REF_KG,
    S=EARTH_HCNS_KRIJT23['S'],
    citation='H/C/N/S: Krijt et al. (2023) PPVII Tables 1+2; O: derived at Sossi 2020 IW+3.5',
    notes=(
        'O is volatile O (atmospheric + dissolved in H2O / CO2 / SO2 / '
        'O2 only), NOT the Krijt et al. Table 2 redox-active O. The '
        'two are inequivalent; see module docstring.'
    ),
)


def derive_earth_volatile_O(T_magma: float = 2000.0, dIW: float = 3.5) -> float:
    """Re-derive Earth's volatile O at given (T_magma, dIW).

    Returns the O_kg_total CALLIOPE's buffered mode reports for the
    Krijt+2023 BSE H/C/N/S target. Useful for confirming the hard-coded
    EARTH_VOLATILE_O_REF_KG has not drifted across CALLIOPE versions.
    """
    import warnings as _warnings

    from calliope.constants import volatile_species
    from calliope.solve import equilibrium_atmosphere

    ddict = {
        'M_mantle': PLANETARY_DEFAULTS['M_mantle'],
        'gravity': PLANETARY_DEFAULTS['gravity'],
        'radius': PLANETARY_DEFAULTS['radius'],
        'Phi_global': 1.0,
        'T_magma': T_magma,
        'fO2_shift_IW': dIW,
    }
    for sp in volatile_species:
        ddict[f'{sp}_included'] = 1
        ddict[f'{sp}_initial_bar'] = 0.0
    with _warnings.catch_warnings():
        _warnings.simplefilter('ignore')
        out = equilibrium_atmosphere(
            EARTH_HCNS_KRIJT23,
            ddict,
            hide_warnings=True,
            print_result=False,
        )
    return float(out['O_kg_total'])


def scale_inventory(base: Inventory, factor: float, name: str | None = None) -> Inventory:
    """Multiply every element budget of `base` by `factor`.

    Useful for parameter sweeps: 0.1x to 10x the canonical Earth budget
    spans the realistic range of volatile-poor to volatile-rich planets.
    """
    if factor <= 0:
        raise ValueError(f'factor must be > 0, got {factor}')
    return Inventory(
        name=name or f'{base.name} x{factor:g}',
        H=base.H * factor,
        C=base.C * factor,
        N=base.N * factor,
        O=base.O * factor,
        S=base.S * factor,
        citation=base.citation + f' (scaled x{factor:g})',
        notes=base.notes,
    )


def scale_O(base: Inventory, O_factor: float, name: str | None = None) -> Inventory:
    """Scale only the O budget. Used by Fig 3 grid: H/C/N/S fixed, O varied."""
    if O_factor <= 0:
        raise ValueError(f'O_factor must be > 0, got {O_factor}')
    return Inventory(
        name=name or f'{base.name} (O x{O_factor:g})',
        H=base.H,
        C=base.C,
        N=base.N,
        O=base.O * O_factor,
        S=base.S,
        citation=base.citation + f' (O scaled x{O_factor:g})',
        notes=base.notes,
    )


PLANETARY_DEFAULTS = {
    'M_mantle': 4.03e24,
    'gravity': 9.81,
    'radius': 6.371e6,
    'Phi_global': 1.0,
    'M_planet': 5.972e24,
    'core_mass_fraction': 0.325,
}
