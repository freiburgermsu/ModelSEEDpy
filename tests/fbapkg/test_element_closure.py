# -*- coding: utf-8 -*-
"""ElementClosurePkg on toy models: a CO2 pump, the biomass loophole, exemptions."""
import cobra
import pytest

from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.fbapkg.elementclosurepkg import formula_atoms


def _toy(pump=True, biomass=True):
    """S (C1) can be taken up at 1. Two reactions form a pump: r1 A -> B + S balances, r2 B -> A
    does not (creates one C), so together they make S from nothing. Biomass consumes a
    placeholder-formula polymer P (formula 'R') made from 2 S."""
    m = cobra.Model("toy")
    M = lambda i, f, q=0: cobra.Metabolite(i, formula=f, charge=q, compartment="c0")
    mets = {"S": M("S_c0", "CO2"), "A": M("A_c0", "C2H4O2"), "B": M("B_c0", "CH4"),
            "P": M("P_c0", "R"), "T": M("cpd11416_c0", "")}
    # make B + S = A in atoms so r1 balances: B = C2H4O2 - CO2 = CH4 (hydrogen differs; give A
    # C2H4O2 and B CH4O0 -> pick A = C2H4O2, B = CH4 means r1 balances C and H but not O:
    # set B = C1H4 and S = CO2, A = C2H4O2 -> r1: C2H4O2 -> CH4 + CO2 : C2 H4 O2 -> C2 H4 O2 ok
    m.add_metabolites(list(mets.values()))

    def rxn(rid, st, lb=0, ub=1000):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[k]: v for k, v in st.items()})
        m.add_reactions([r])

    rxn("EX_S_e0", {"S": -1}, lb=-1, ub=1000)
    if pump:
        rxn("r1", {"A": -1, "B": 1, "S": 1})          # balanced
        rxn("r2", {"B": -1, "A": 1})                  # creates CO2 worth of atoms
    if biomass:
        rxn("Psyn", {"S": -2, "P": 1})                # P's content = 2 CO2
        rxn("bio1", {"P": -1, "T": 1})
        rxn("SK_cpd11416_c0", {"T": -1})
    return m


def _close(m, **params):
    return MSPackageManager.get_pkg_mgr(m).getpkg("ElementClosurePkg").build_package(params)


def _max(m, rid):
    m.objective = rid
    v = m.slim_optimize()
    return 0.0 if v != v else v


def test_formula_parsing():
    assert formula_atoms("C4H5N2O3R2") is None
    assert formula_atoms("C4H5N2O3R2", strip_r=True) == {"C": 4, "H": 5, "N": 2, "O": 3}
    assert formula_atoms("C62H90CoN13O14PR", strip_r=True)["Co"] == 1
    assert formula_atoms("RbCl") == {"Rb": 1, "Cl": 1}          # Rb is not an R group
    assert formula_atoms("") is None


def test_pump_exports_carbon_without_closure_and_not_with_it():
    m = _toy(biomass=False)
    assert _max(m, "EX_S_e0") == pytest.approx(1000)            # carbon out of nothing
    m = _toy(biomass=False)
    _close(m)
    assert _max(m, "EX_S_e0") == pytest.approx(0, abs=1e-7)


def test_biomass_loophole_is_closed():
    # without closure the pump feeds biomass without limit
    assert _max(_toy(), "bio1") == pytest.approx(500.5)        # (1000 pump + 1 uptake) / 2
    # exchange-only closure (biomass not counted) still lets pump carbon into biomass
    m = _toy()
    _close(m, biomass_reactions=[])
    assert _max(m, "bio1") == pytest.approx(500.5)
    # with biomass counted (P resolved to 2 CO2 through Psyn): one S in supports 0.5 biomass
    m = _toy()
    info = _close(m)
    assert info["pseudo_compounds"]["P_c0"] == "synthesis Psyn"
    assert _max(m, "bio1") == pytest.approx(0.5)


def test_real_pathway_is_untouched():
    m = _toy(pump=False)
    _close(m)
    assert _max(m, "bio1") == pytest.approx(0.5)                # same as with no closure
    m2 = _toy(pump=False)
    assert _max(m2, "bio1") == pytest.approx(0.5)


def test_exempt_reaction_may_create():
    m = _toy(biomass=False)
    _close(m, exempt_reactions={"r2": "known error, pending curation"})
    assert _max(m, "EX_S_e0") == pytest.approx(1000)


def test_creation_report_and_curation_targets():
    m = _toy(biomass=False)
    pkg = MSPackageManager.get_pkg_mgr(m).getpkg("ElementClosurePkg")
    pkg.build_package({"elements": ("C", "O")})
    m.remove_cons_vars([pkg.constraints["closure"][k] for k in list(pkg.constraints["closure"])])
    m.objective = "EX_S_e0"
    sol = m.optimize()
    made = pkg.creation(sol.fluxes)
    assert made["C"] == pytest.approx(1000) and made["O"] == pytest.approx(2000)
    assert set(pkg.imbalanced_reactions()) == {"r2"}


def test_all_elements_and_charge_rows():
    m = _toy()
    info = _close(m, elements="all", charge=True)
    assert set(info["rows"]) == {"C", "O", "charge"}
