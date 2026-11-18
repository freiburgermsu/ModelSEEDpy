# -*- coding: utf-8 -*-
from modelseedpy.core.msgapfill import MSGapfill
from tests.core.test_msgapfill import (  # noqa: F401
    get_model,
    media_glucose_aerobic,
    template,
)


def test_binary_check_reads_the_minimised_solution(
    template, get_model, media_glucose_aerobic, monkeypatch
):
    """The binary check minimises the number of gapfilled reactions and must then
    read the solution from that minimised optimum, not from the fluxes captured
    before it."""
    model = get_model(["CYTBD_c0", "NADH16_c0", "GLCpts_c0", "O2t_c0"])
    gap_fill = MSGapfill(model, [template])
    gfpkg = gap_fill.gfpkgmgr.getpkg("GapfillingPkg")
    flux_args = []
    original = gfpkg.compute_gapfilled_solution

    def recording(flux_values=None):
        flux_args.append(flux_values)
        return original(flux_values)

    monkeypatch.setattr(gfpkg, "compute_gapfilled_solution", recording)
    result = gap_fill.run_gapfilling(
        media_glucose_aerobic,
        "BIOMASS_Ecoli_core_w_GAM_c0",
        minimum_obj=0.8,
        binary_check=True,
    )
    assert result
    # the last read is the one after the binary minimisation: fresh fluxes
    assert flux_args[-1] is None
    for rxn in ["CYTBD_c0", "NADH16_c0", "GLCpts_c0", "O2t_c0"]:
        assert result["new"].get(rxn) == ">"
