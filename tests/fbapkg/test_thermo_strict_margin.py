# -*- coding: utf-8 -*-
"""The loop-law margin in SimpleThermoPkg, and the units of FullThermoPkg's database errors.

Toy models only -- no database or solver licence beyond the default MILP solver is needed.
"""
import cobra
import pytest

from modelseedpy.fbapkg.mspackagemanager import MSPackageManager


def _model(with_loop=True, with_pump=False):
    """A -> B -> C plus exchanges; optionally an internal cycle B <-> C <-> D <-> B, or a
    two-reaction 'pump' P + Q -> X + R, R -> P + Q... that nets nothing -> X."""
    m = cobra.Model("toy")
    mets = {k: cobra.Metabolite(f"{k}_c0", compartment="c0") for k in "ABCDXPQR"}
    m.add_metabolites(list(mets.values()))

    def rxn(rid, stoich, lb=-1000, ub=1000):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[k]: v for k, v in stoich.items()})
        m.add_reactions([r])
        return r

    rxn("EX_A_e0", {"A": -1}, lb=-10, ub=0)            # A in
    rxn("EX_C_e0", {"C": -1}, lb=0, ub=1000)            # C out
    rxn("R1", {"A": -1, "B": 1})
    rxn("R2", {"B": -1, "C": 1})
    if with_loop:                                       # a cycle with no net conversion
        rxn("L1", {"C": -1, "D": 1})
        rxn("L2", {"D": -1, "B": 1})
    if with_pump:                                       # nothing -> X, through P/Q/R
        rxn("P1", {"P": -1, "Q": -1, "X": 1, "R": 1})
        rxn("P2", {"R": -1, "P": 1, "Q": 1})
        rxn("EX_X_e0", {"X": -1}, lb=0, ub=1000)
    return m


def _max(m, rid, **params):
    if params is not None:
        MSPackageManager.get_pkg_mgr(m).getpkg("SimpleThermoPkg").build_package(params)
    m.objective = rid
    return m.slim_optimize()


@pytest.mark.parametrize("coupling", ["revbin", "tmfa"])
def test_margin_blocks_internal_cycle(coupling):
    plain = _model()
    plain.objective = "L1"
    assert plain.slim_optimize() == pytest.approx(1000)          # the cycle runs freely
    # eps = 0 is the old behaviour: all-zero potentials satisfy every row
    assert _max(_model(), "L1", dg_margin=0.0, coupling=coupling) == pytest.approx(1000)
    # any eps > 0 forbids the cycle; L1 may only carry flux the network actually needs (none)
    assert _max(_model(), "L1", dg_margin=1.0, coupling=coupling) == pytest.approx(0, abs=1e-6)


@pytest.mark.parametrize("coupling", ["revbin", "tmfa"])
def test_margin_keeps_real_pathways(coupling):
    # A -> B -> C -> out is downhill all the way; the margin must not cost it any flux
    assert _max(_model(with_loop=False), "EX_C_e0", dg_margin=1.0, coupling=coupling) \
        == pytest.approx(10)


@pytest.mark.parametrize("coupling", ["revbin", "tmfa"])
def test_margin_cannot_stop_mass_creation(coupling):
    # P1 + P2 net: nothing -> X. It is not a cycle (X leaves), so no loop law excludes it.
    # Documented limitation: this needs mass balance, not thermodynamics.
    m = _model(with_loop=False, with_pump=True)
    assert _max(m, "EX_X_e0", dg_margin=1.0, coupling=coupling) == pytest.approx(1000)


def test_tmfa_leaves_idle_reactions_unconstrained():
    m = _model()
    pkg = MSPackageManager.get_pkg_mgr(m).getpkg("SimpleThermoPkg")
    pkg.build_package({"dg_margin": 1.0, "coupling": "tmfa"})
    # pin R1's dgp at exactly 0: under tmfa R1 can then carry no flux, but the model stays
    # feasible; under revbin the same pin is infeasible because every reaction needs |dgp|>=eps
    pkg.variables["dgp"]["R1"].lb = pkg.variables["dgp"]["R1"].ub = 0
    m.objective = "EX_C_e0"
    assert m.slim_optimize() == pytest.approx(0, abs=1e-6)
    r = _model()
    rp = MSPackageManager.get_pkg_mgr(r).getpkg("SimpleThermoPkg")
    rp.build_package({"dg_margin": 1.0, "coupling": "revbin"})
    rp.variables["dgp"]["R1"].lb = rp.variables["dgp"]["R1"].ub = 0
    r.objective = "EX_C_e0"
    r.slim_optimize()                    # nan, not an exception
    assert r.solver.status == "infeasible"


def test_build_package_without_arguments():
    # MetaboFBAPkg calls build_package() with no parameters
    m = _model()
    MSPackageManager.get_pkg_mgr(m).getpkg("SimpleThermoPkg").build_package()
    m.objective = "L1"
    assert m.slim_optimize() == pytest.approx(0, abs=1e-6)


def test_loader_skips_equilibrator_no_estimate_marker():
    from modelseedpy.biochem.modelseed_biochem import _thermodynamics
    sam = {"thermodynamics": {"eQuilibrator": [132.52, 23901.43],
                              "dGPredictor": [67.83, 8.41, 0.76]}}
    assert _thermodynamics(sam) == (67.83, 8.41, "dGPredictor")
    co2 = {"thermodynamics": {"eQuilibrator": [-92.26, 0.76], "dGPredictor": [-91.91, 2.09]}}
    assert _thermodynamics(co2) == (-92.26, 0.76, "eQuilibrator")
    unknown = {"thermodynamics": {"eQuilibrator": [27.8, 12482.4],
                                  "Group contribution": [1e7, 1e7]}}
    assert _thermodynamics(unknown) == (None, None, None)


def test_database_error_is_converted_to_kj(monkeypatch):
    from modelseedpy.fbapkg import fullthermopkg as ft

    class Cpd:
        def __init__(self, err):
            self.delta_g_error = err

    class Compounds(dict):
        def get_by_id(self, k):
            return self[k]

    class DB:
        compounds = Compounds(cpd1=Cpd(2.0), cpd2=Cpd(23901.43), cpd3=Cpd(None))

    monkeypatch.setattr("modelseedpy.biochem.modelseed_biochem.ModelSEEDBiochem.get",
                        staticmethod(lambda path=None: DB()))
    pkg = ft.FullThermoPkg.__new__(ft.FullThermoPkg)
    pkg.parameters = {"modelseed_db_path": None}
    assert pkg.database_error("cpd1") == pytest.approx(2.0 * 4.184)   # kcal -> kJ
    assert pkg.database_error("cpd2") is None                         # no-estimate marker
    assert pkg.database_error("cpd3") is None
    assert pkg.database_error("cpdX") is None
