# -*- coding: utf-8 -*-
from cobra import Metabolite, Model, Reaction

from modelseedpy.core.msmodelutl import MSModelUtil


def test_exchange_list_matches_exchange_prefixes_only():
    model = Model("exchanges")
    met = Metabolite("cpd00001_e0", compartment="e0")
    for rxn_id in ["EX_cpd00001_e0", "EXF_cpd00001_e0", "SK_EX_cpd00001_e0"]:
        rxn = Reaction(rxn_id)
        rxn.add_metabolites({met: -1})
        model.add_reactions([rxn])
    ids = sorted(rxn.id for rxn in MSModelUtil(model).exchange_list())
    assert ids == ["EXF_cpd00001_e0", "EX_cpd00001_e0"]
