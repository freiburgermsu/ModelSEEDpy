# -*- coding: utf-8 -*-
from modelseedpy.core.mstemplate import MSTemplateMetabolite, NewModelTemplateRole


def test_template_metabolite_html_shows_values():
    html = MSTemplateMetabolite("cpd00001", "H2O", "Water")._repr_html_()
    assert "cpd00001" in html and "Water" in html
    assert "{self.id}" not in html


def test_template_role_html_shows_values():
    html = NewModelTemplateRole("ftr01", "Glucokinase")._repr_html_()
    assert "ftr01" in html and "Glucokinase" in html
    assert "{self.id}" not in html
