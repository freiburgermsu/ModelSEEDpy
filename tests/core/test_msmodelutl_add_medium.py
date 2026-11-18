# -*- coding: utf-8 -*-
from modelseedpy.core.msmedia import MSMedia
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.test_data.mock_data import mock_model_ecoli_core


def _model_with_seed_exchanges():
    model = mock_model_ecoli_core(seed=True)
    model.reactions.get_by_id("EX_glc__D_e").id = "EX_cpd00027_e0"
    model.reactions.get_by_id("EX_o2_e").id = "EX_cpd00007_e0"
    model.repair()
    return model


def test_add_medium_keeps_the_media_bounds():
    model = _model_with_seed_exchanges()
    media = MSMedia.from_dict({"cpd00027": 10, "cpd00007": 20})
    medium = MSModelUtil(model).add_medium(media)
    assert medium["EX_cpd00027_e0"] == 10
    assert medium["EX_cpd00007_e0"] == 20


def test_add_medium_uniform_uptake_still_overrides():
    model = _model_with_seed_exchanges()
    media = MSMedia.from_dict({"cpd00027": 10, "cpd00007": 20})
    medium = MSModelUtil(model).add_medium(media, uniform_uptake=5)
    assert medium["EX_cpd00027_e0"] == 5
    assert medium["EX_cpd00007_e0"] == 5
