# -*- coding: utf-8 -*-
"""MSModelUtil.add_ms_reaction against a small stand-in for the ModelSEED database.

The stand-in's reactions are built by modelseed_biochem.build_modelseed_reaction, the
token construction the database loader uses, so their metabolites carry compartment
index ids such as cpd00001_0, as in the real database. The model is the mock
e_coli_core with ModelSEED compound ids (cpd00001_c0, ...).
"""
from types import SimpleNamespace

import pytest
from cobra import DictList

from modelseedpy.biochem.modelseed_biochem import build_modelseed_reaction
from modelseedpy.biochem.modelseed_compound import ModelSEEDCompound2
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.test_data.mock_data import mock_model_ecoli_core

# an id the mock model lacks: it already holds rxn00062_c0, which cobra would ignore
RXN_ID = "rxn99999"
# (compound, compartment index, coefficient): ATP + H2O -> ADP + Pi + H+, as rxn00062
ATP_HYDROLYSIS = [
    ("cpd00002", 0, -1),
    ("cpd00001", 0, -1),
    ("cpd00008", 0, 1),
    ("cpd00009", 0, 1),
    ("cpd00067", 0, 1),
]


def _database(stoichiometry):
    compounds = {
        cpd: ModelSEEDCompound2(cpd, name=cpd.upper()) for cpd, _, _ in stoichiometry
    }
    reaction = build_modelseed_reaction(
        {
            "id": RXN_ID,
            "name": "test reaction",
            "reversibility": ">",
            "stoichiometry": [
                {"compound": cpd, "compartment": index, "coefficient": coef}
                for cpd, index, coef in stoichiometry
            ],
        },
        {},
        {},
        {},
        {},
        compounds,
    )
    return SimpleNamespace(reactions=DictList([reaction]))


def _refuse_from_local(path):
    raise AssertionError(f"from_local({path!r}) was called")


def test_default_database_is_modelseedbiochem_get(monkeypatch):
    db = _database(ATP_HYDROLYSIS)
    monkeypatch.setattr(
        "modelseedpy.biochem.modelseed_biochem.ModelSEEDBiochem.get",
        staticmethod(lambda *args, **kwargs: db),
    )
    monkeypatch.setattr("modelseedpy.biochem.from_local", _refuse_from_local)
    model = mock_model_ecoli_core(seed=True)
    (rxn,) = MSModelUtil(model).add_ms_reaction({RXN_ID: "c0"})
    assert rxn is model.reactions.get_by_id(RXN_ID + "_c0")


def test_database_path_is_loaded_with_from_local(monkeypatch):
    db = _database(ATP_HYDROLYSIS)
    paths = []
    monkeypatch.setattr(
        "modelseedpy.biochem.from_local", lambda path: paths.append(path) or db
    )
    model = mock_model_ecoli_core(seed=True)
    # msdb_path stays the second positional parameter
    MSModelUtil(model).add_ms_reaction({RXN_ID: "c0"}, "/path/to/ModelSEEDDatabase")
    assert paths == ["/path/to/ModelSEEDDatabase"]
    assert RXN_ID + "_c0" in model.reactions


@pytest.mark.parametrize("keyword", ["compartment_trans", "comp_trans", None])
def test_compartment_translation(keyword):
    # Pi from compartment index 0 to index 1, translated to c0 and p0
    db = _database([("cpd00009", 0, -1), ("cpd00009", 1, 1)])
    model = mock_model_ecoli_core(seed=True)
    mdlutl = MSModelUtil(model)
    if keyword:
        (rxn,) = mdlutl.add_ms_reaction(
            {RXN_ID: "c0"}, msdb_object=db, **{keyword: ["c0", "p0"]}
        )
    else:  # positional: rxn_dict, msdb_path, msdb_object, comp_trans
        (rxn,) = mdlutl.add_ms_reaction({RXN_ID: "c0"}, None, db, ["c0", "p0"])
    assert {met.id: coef for met, coef in rxn.metabolites.items()} == {
        "cpd00009_c0": -1,
        "cpd00009_p0": 1,
    }
    assert model.metabolites.get_by_id("cpd00009_p0").compartment == "p0"


def test_reuses_model_metabolites_by_compound_id():
    model = mock_model_ecoli_core(seed=True)
    before = len(model.metabolites)
    (rxn,) = MSModelUtil(model).add_ms_reaction(
        {RXN_ID: "c0"}, msdb_object=_database(ATP_HYDROLYSIS)
    )
    assert rxn is model.reactions.get_by_id(RXN_ID + "_c0")
    assert {met.id: coef for met, coef in rxn.metabolites.items()} == {
        f"{cpd}_c0": coef for cpd, _, coef in ATP_HYDROLYSIS
    }
    # the reaction is wired to the model's own metabolites, so nothing was added
    for met in rxn.metabolites:
        assert met is model.metabolites.get_by_id(met.id)
    assert len(model.metabolites) == before


def test_new_compound_is_named_from_compound_id():
    model = mock_model_ecoli_core(seed=True)
    before = len(model.metabolites)
    (rxn,) = MSModelUtil(model).add_ms_reaction(
        {RXN_ID: "c0"}, msdb_object=_database([("cpd00001", 0, -1), ("cpd90000", 0, 1)])
    )
    assert sorted(met.id for met in rxn.metabolites) == ["cpd00001_c0", "cpd90000_c0"]
    assert model.metabolites.get_by_id("cpd90000_c0").compartment == "c0"
    assert len(model.metabolites) == before + 1
