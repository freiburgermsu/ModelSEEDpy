# -*- coding: utf-8 -*-
"""ElementClosurePkg -- the network may not create chemical elements.

For each element E one linear row, over the reactions that cross the system boundary:

    sum over boundary and biomass reactions r of  -imb_E(r) * v_r   <=  0

where imb_E(r) is the net number of E atoms reaction r produces per unit flux. An exchange
that exports a compound "destroys" its atoms (imb < 0), so the row reads: E leaving through
exchanges, sinks and demands, plus E retained in biomass, may not exceed E entering.

Why the boundary is enough. At steady state S v = 0, so for ANY assignment a of atom counts to
metabolites, a^T S v = 0: the net creation of E by internal reactions equals minus its net
creation by boundary and biomass reactions. Forbidding net creation at the boundary therefore
forbids it inside, whatever the internal mechanism -- including reactions whose own balance
cannot even be computed because a reactant carries an R group or a placeholder formula. One row
per element, no binaries.

What it catches and what it does not.
  - Mass-imbalanced reactions or combinations that create atoms -- e.g. two versions of one
    enzymatic step that differ by a CO2, which together emit CO2 from nothing. Thermodynamic
    (loop-law) constraints cannot catch these: such a pathway is not a cycle, and creating CO2
    "from nothing" looks downhill to every formation-energy model.
  - NOT charge-only imbalances that create electrons with no atoms (e.g. carrier pairs that
    differ in how many ferredoxins they reduce). Set charge=True to add a charge row; expect it
    to expose every charge-imbalanced reaction a model relies on.
  - It forbids net creation, not every imbalance: a reaction that creates an element may run if
    another destroys as much (the net is what reaches the boundary). imbalanced_reactions()
    lists the culprits so they can be curated.

Biomass. Biomass is a sink that keeps atoms. Its content is computed from its components; the
biomass token (cpd11416, which every biomass reaction emits and sinks remove) carries none, so
nothing is counted twice and member-level biomass sinks in a community are covered. Macromolecule
pseudo-compounds -- Protein, DNA, RNA, ACP -- have placeholder formulas (C4H5N2O3R2 for Protein),
so each is given the content of the reaction that synthesises it from formulated precursors. A
carrier with no producer (apo-ACP, which biomass releases as it consumes ACP) is zero, which makes
the pair count only the moiety actually incorporated.

Atom assignment (one, used everywhere): biomass token -> 0; formula without R -> its atoms;
biomass pseudo-compound -> content of its synthesis reaction; any other formula with R -> its
explicit atoms (R taken as the same unspecified moiety on both sides of a reaction); no formula
-> 0, and listed in `unformulated`.

Parameters (build_package):
    elements          iterable of element symbols, or "all" (every element the boundary and
                      biomass compounds contain). Default C, H, O, N, P, S.
    charge            also close charge. Default False.
    sense             "le" (default): no net creation, destruction allowed -- biomass may keep
                      atoms, and so may sinks. "eq": neither.
    exempt_reactions  {reaction id: reason} -- reactions whose own imbalance is allowed (e.g. a
                      known stoichiometric error awaiting curation that growth depends on). Every
                      other reaction stays closed. Their metabolites must have formulas.
    biomass_reactions ids; default: every reaction that produces a biomass token.
    biomass_tokens    ModelSEED ids of biomass tokens. Default ("cpd11416",).
"""
from __future__ import absolute_import

import logging
import re

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

logger = logging.getLogger(__name__)

DEFAULT_ELEMENTS = ("C", "H", "O", "N", "P", "S")
_ATOM = re.compile(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)?")
_RGROUP = re.compile(r"R(?![a-z])(\d+(?:\.\d+)?)?")      # R, not Rb/Re/Rh/Ru/Rn/Ra/Rf/Rg


def formula_atoms(formula, strip_r=False):
    """Atom counts of a formula. None if there is none, or if it carries an R group and
    strip_r is False; with strip_r the R groups are dropped and the explicit atoms counted."""
    if not formula or formula in ("null", "noformula"):
        return None
    if _RGROUP.search(formula):
        if not strip_r:
            return None
        formula = _RGROUP.sub("", formula)
    out = {}
    for el, n in _ATOM.findall(formula):
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def _base_id(met):
    return met.id.rsplit("_", 1)[0] if "_" in met.id else met.id


class ElementClosurePkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "element closure", {}, {"closure": "string"})
        self.content = {}            # met id -> (atoms dict, charge)
        self.how = {}                # met id -> how its content was assigned
        self.biomass = []
        self.boundary = []
        self.terms = {}              # element -> {reaction id: -imb (coefficient on net flux)}

    # ------------------------------------------------------------------ atom assignment ------
    def _parsed(self, met):
        atoms = formula_atoms(met.formula)
        if atoms is not None:
            return atoms, float(met.charge or 0)
        return None

    def build_contents(self):
        p = self.parameters
        tokens = set(p["biomass_tokens"])
        for met in self.model.metabolites:
            if _base_id(met) in tokens:
                self.content[met.id], self.how[met.id] = ({}, 0.0), "biomass token"
        if p["biomass_reactions"] is None:
            self.biomass = [r for r in self.model.reactions
                            if any(_base_id(m) in tokens for m in r.products)]
        else:
            self.biomass = [self.model.reactions.get_by_id(r) for r in p["biomass_reactions"]]
        bio_ids = {r.id for r in self.biomass}

        def known(met):
            if met.id in self.content:
                return self.content[met.id]
            v = self._parsed(met)
            if v is not None:
                self.content[met.id], self.how[met.id] = v, "formula"
            return v

        # biomass components with placeholder formulas: give each the content of the reaction
        # that synthesises it; carriers with no producer are zero; repeat until nothing moves
        unknown = {m for b in self.biomass for m in b.metabolites if known(m) is None}
        for _ in range(50):
            moved = False
            for u in sorted(unknown, key=lambda m: m.id):
                for r in u.reactions:
                    if r.id in bio_ids or r.boundary or r.metabolites[u] <= 0:
                        continue
                    others = [m for m in r.metabolites if m is not u]
                    vals = [known(m) for m in others]
                    if any(v is None for v in vals):
                        continue
                    s = r.metabolites[u]
                    atoms, q = {}, 0.0
                    for m, (a, c) in zip(others, vals):
                        for el, n in a.items():
                            atoms[el] = atoms.get(el, 0.0) + r.metabolites[m] * n
                        q += r.metabolites[m] * c
                    self.content[u.id] = ({el: -n / s for el, n in atoms.items()
                                           if abs(n) > 1e-12}, -q / s)
                    self.how[u.id] = f"synthesis {r.id}"
                    unknown.discard(u)
                    moved = True
                    break
            if moved:
                continue
            dead = [u for u in unknown if not any(
                r.metabolites[u] > 0 and r.id not in bio_ids and not r.boundary
                for r in u.reactions)]
            for u in (dead or list(unknown)):
                self.content[u.id] = ({}, 0.0)
                self.how[u.id] = "zero: carrier with no producer" if dead else \
                    "zero: unresolved biomass component"
                unknown.discard(u)
                if not dead:
                    logger.warning(f"ElementClosurePkg: biomass component {u.id} has no "
                                   f"resolvable content; counted as zero")
            if not unknown:
                break

    def atoms_of(self, met):
        """Content of a metabolite under this package's single assignment."""
        if met.id in self.content:
            return self.content[met.id]
        v = self._parsed(met)
        if v is not None:
            how = "formula"
        else:
            a = formula_atoms(met.formula, strip_r=True)
            v, how = ((a, float(met.charge or 0)), "explicit atoms (R dropped)") \
                if a is not None else (({}, float(met.charge or 0)), "zero: no formula")
        self.content[met.id], self.how[met.id] = v, how
        return v

    def imbalance(self, rxn, key):
        """Net atoms of `key` (an element, or 'charge') the reaction produces per unit flux."""
        tot = 0.0
        for met, c in rxn.metabolites.items():
            atoms, q = self.atoms_of(met)
            tot += c * (q if key == "charge" else atoms.get(key, 0.0))
        return tot

    # ------------------------------------------------------------------ constraints ----------
    def build_package(self, parameters=None):
        terms = self.closure_terms(parameters)
        p = self.parameters
        for key, row in terms.items():
            coef = {}
            for rid, c in row.items():
                r = self.model.reactions.get_by_id(rid)
                coef[r.forward_variable] = coef.get(r.forward_variable, 0.0) + c
                coef[r.reverse_variable] = coef.get(r.reverse_variable, 0.0) - c
            self.build_constraint("closure", 0 if p["sense"] == "eq" else None, 0, coef, key)
        return self.summary()

    def closure_terms(self, parameters=None):
        """The closure rows without building them: {element: {reaction id: coefficient}}, each
        row reading sum(coefficient * net flux) <= 0 (== 0 with sense "eq").

        For solvers driven outside optlang -- e.g. a model too large to load into an optlang
        interface, solved from its stoichiometric matrix -- add these rows directly.
        build_package builds exactly these.
        """
        self.validate_parameters(
            parameters or {},
            [],
            {
                "elements": DEFAULT_ELEMENTS,
                "charge": False,
                "sense": "le",
                "exempt_reactions": {},
                "biomass_reactions": None,
                "biomass_tokens": ("cpd11416",),
            },
        )
        p = self.parameters
        if p["sense"] not in ("le", "eq"):
            raise ValueError("sense must be 'le' or 'eq'")
        exempt = p["exempt_reactions"]
        if not isinstance(exempt, dict):
            exempt = {r: "" for r in exempt}
        self.build_contents()
        bio_ids = {r.id for r in self.biomass}
        self.boundary = [r for r in self.model.reactions
                         if r.id not in bio_ids and (r.boundary or
                                                     r.id.startswith(("EX_", "SK_", "DM_")))]
        rows = list(self.boundary) + list(self.biomass)
        elements = p["elements"]
        if elements == "all":
            found = set()
            for r in rows:
                for met in r.metabolites:
                    found.update(self.atoms_of(met)[0])
            elements = sorted(found)
        keys = list(elements) + (["charge"] if p["charge"] else [])
        exempt_rxns = []
        for rid, why in exempt.items():
            if rid not in self.model.reactions:
                continue
            r = self.model.reactions.get_by_id(rid)
            bad = [m.id for m in r.metabolites if self._parsed(m) is None
                   and m.id not in self.content]
            if bad:
                raise ValueError(f"exempt reaction {rid} has metabolites without a formula: {bad}")
            exempt_rxns.append(r)
            logger.info(f"ElementClosurePkg: {rid} exempt ({why})")
        out = {}
        for key in keys:
            terms = {}
            for r in rows:
                b = self.imbalance(r, key)
                if abs(b) > 1e-12:
                    terms[r.id] = -b
            for r in exempt_rxns:             # allow this reaction's own creation of `key`
                b = self.imbalance(r, key)
                if abs(b) > 1e-12:
                    terms[r.id] = terms.get(r.id, 0.0) - b
            self.terms[key] = terms
            out[key] = dict(terms)
        return out

    # ------------------------------------------------------------------ reporting ------------
    def summary(self):
        how = {}
        for v in self.how.values():
            k = v.split(" ")[0] if v.startswith("synthesis") else v
            how[k] = how.get(k, 0) + 1
        return {"rows": sorted(self.terms), "boundary_reactions": len(self.boundary),
                "biomass_reactions": len(self.biomass),
                "pseudo_compounds": {m: h for m, h in self.how.items()
                                     if h.startswith(("synthesis", "zero"))},
                "assignment_counts": how}

    def creation(self, fluxes):
        """Net creation of each closed quantity by internal reactions outside the exemptions,
        from a flux mapping (reaction id -> flux). This is the left-hand side of the row:
        positive means the solution creates it, and the row keeps it <= 0. By a^T S v = 0 it
        equals the sum over internal reactions of imbalance x flux."""
        out = {}
        for key, terms in self.terms.items():
            out[key] = sum(c * float(fluxes[rid]) for rid, c in terms.items() if rid in fluxes)
        return out

    def imbalanced_reactions(self, keys=None, tol=1e-9):
        """Internal reactions that create or destroy a closed quantity, under formulas alone
        (reactions involving an R-group or unformulated metabolite are left out: their balance
        is not defined). These are the curation targets."""
        keys = keys or list(self.terms)
        bio_ids = {r.id for r in self.biomass}
        bnd = {r.id for r in self.boundary}
        out = {}
        for r in self.model.reactions:
            if r.id in bio_ids or r.id in bnd:
                continue
            if any(self._parsed(m) is None and not str(self.how.get(m.id, "")).startswith(
                    ("synthesis", "biomass token")) for m in r.metabolites):
                continue
            imb = {k: self.imbalance(r, k) for k in keys}
            imb = {k: v for k, v in imb.items() if abs(v) > tol}
            if imb:
                out[r.id] = imb
        return out
