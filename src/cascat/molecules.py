"""Small, extensible catalogue; names encode bonding order, not just formula."""

# Public name: (ACAT binding name, reaction family).
CATALOGUE = {
    "H": ("H", "HER/CO2RR/OER"),
    "H2": ("H2", "HER"),
    "H2O": ("OH2", "HER/CO2RR/OER"),
    "O": ("O", "OER"),
    "OH": ("OH", "HER/CO2RR/OER"),
    "OOH": ("OOH", "OER"),
    "O2": ("O2", "OER"),
    "H2O2": ("OHOH", "OER"),
    "CO2": ("CO2", "CO2RR"),
    "COOH": ("COOH", "CO2RR"),
    "HCOO": ("CHOO", "CO2RR"),
    "HCOOH": ("CHOOH", "CO2RR"),
    "CO": ("CO", "CO2RR"),
    "CHO": ("CHO", "CO2RR"),
    "COH": ("COH", "CO2RR"),
    "CHOH": ("CHOH", "CO2RR"),
    "CH2O": ("CH2O", "CO2RR"),
    "CH2OH": ("CH2OH", "CO2RR"),
    "CH3O": ("OCH3", "CO2RR"),
    "CH3OH": ("CH3OH", "CO2RR"),
    "C": ("C", "CO2RR"),
    "CH": ("CH", "CO2RR"),
    "CH2": ("CH2", "CO2RR"),
    "CH3": ("CH3", "CO2RR"),
    "CH4": ("CH4", "CO2RR"),
    "CH2CO": ("CH2CO", "CO2RR/C2"),
    "CH3CO": ("CH3CO", "CO2RR/C2"),
    "CH3COOH": ("CH3COOH", "CO2RR/C2"),
}


def configure(engine, public_name, minimum_distance=None):
    """Register supported ACAT templates without modifying installed packages."""
    from collections import Counter
    from acat import settings

    name = CATALOGUE[public_name][0]
    template = settings.adsorbate_molecule(name)
    if template is None:
        raise ValueError(f"No geometry template for {public_name}")
    if name not in settings.adsorbate_list:
        settings.adsorbate_list.append(name)
        settings.monodentate_adsorbate_list.append(name)
        settings.adsorbate_formulas[name] = "".join(template.get_chemical_symbols())
    engine.ADSORBATE_STOICHIOMETRY[name] = dict(Counter(template.get_chemical_symbols()))
    engine.ADSORBATE_ANCHOR[name] = template[0].symbol
    engine.MIN_ADSORBATE_DISTANCE[name] = (
        minimum_distance if minimum_distance is not None
        else (1.8 if name == "H" else 2.2)
    )
    engine.CLEAR_DIMER_DISTANCE.setdefault(name, None)
    return name
