# -*- coding: utf-8 -*-
import logging
import os

from modelseedpy.core import msmodelutl
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.test_data.mock_data import mock_model_ecoli_core


class _Media:
    id = "infeasible_media"


def _infeasible_condition_result():
    model = mock_model_ecoli_core(seed=True)
    model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = 1000
    condition = {"media": _Media(), "threshold": 0.1, "is_max_threshold": True}
    return MSModelUtil(model).test_single_condition(
        condition, apply_condition=False, model=model
    )


def test_infeasible_condition_writes_no_file_by_default(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert _infeasible_condition_result() is False
    assert os.listdir(tmp_path) == []


def test_infeasible_condition_lp_is_written_when_debugging(
    tmp_path, monkeypatch, caplog
):
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.DEBUG, logger=msmodelutl.logger.name)
    assert _infeasible_condition_result() is False
    assert os.listdir(tmp_path) == ["infeasible_media-Testing-Infeasible.lp"]
