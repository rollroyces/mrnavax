"""LNP composition advisor.

Recommends an LNP formulation based on:
- target tissue (lung, liver, spleen, muscle, tumor)
- cargo type (mRNA, saRNA, circRNA, siRNA, Cas9 mRNA)
- therapeutic intent (cancer vaccine, gene editing, protein replacement)

The recommendation is a **rule-based shortlist** of published or ML-discovered
ionizable lipids + helper lipid ratios. The ML layer in production would be
something like the directed message-passing neural network in
Witten et al. *Nat Biotech* 43, 1790–1799 (2025) — the LION nanoparticle library,
or Li et al. *Nat Mater* 23, 1002 (2024) combinatorial-chemistry accelerated
ionizable-lipid discovery. This module is the **interface layer** above that
model: humans still pick from a shortlist.

References
----------
Witten et al. 2025 (AI LNP for lung):    https://www.nature.com/articles/s41587-024-02490-y
Li et al. 2024 (combinatorial+ML):        https://www.nature.com/articles/s41563-024-01833-5
Hou et al. 2021 (LNP review):             Nat Rev Mater 6, 1078
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

# ---------- curated table of published LNP formulations ---------------------


@dataclass
class LNPPreset:
    name: str
    ionizable_lipid: str
    helper_lipid: str
    cholesterol_pct: float
    peg_lipid: str
    peg_mol_pct: float
    ionizable_mol_pct: float
    helper_mol_pct: float
    n_p_ratio: float  # N/P ratio
    source: str
    notes: str

    def as_dict(self) -> dict:
        return asdict(self)


PRESETS: dict[str, LNPPreset] = {
    # Moderna's benchmark
    "SM-102": LNPPreset(
        name="SM-102 (Moderna)",
        ionizable_lipid="SM-102",
        helper_lipid="DSPC",
        cholesterol_pct=38.5,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.5,
        ionizable_mol_pct=50.0,
        helper_mol_pct=10.0,
        n_p_ratio=6.0,
        source="Hou et al. 2021 / Spikevax FDA label",
        notes="Clinical benchmark for mRNA vaccines; hepatic + splenic delivery.",
    ),
    "ALC-0315": LNPPreset(
        name="ALC-0315 (Pfizer/BioNTech)",
        ionizable_lipid="ALC-0315",
        helper_lipid="DSPC",
        cholesterol_pct=42.7,
        peg_lipid="ALC-0159",
        peg_mol_pct=1.5,
        ionizable_mol_pct=46.3,
        helper_mol_pct=9.4,
        n_p_ratio=6.0,
        source="Hou et al. 2021 / Comirnaty EMA EPAR",
        notes="Clinical benchmark; broader tropism than SM-102 in some tissues.",
    ),
    # C12-200 — classic hepatocyte-delivery ionizable
    "C12-200": LNPPreset(
        name="C12-200 (hepatocyte)",
        ionizable_lipid="C12-200",
        helper_lipid="DOPE",
        cholesterol_pct=28.0,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.5,
        ionizable_mol_pct=50.0,
        helper_mol_pct=10.0,
        n_p_ratio=4.5,
        source="Love et al. PNAS 2010 / numerous reviews",
        notes="Strong hepatocyte uptake; widely used for hepatic gene editing.",
    ),
    # FO-32 / FO-35 — Witten et al. 2025 lung-delivery lipids
    "FO-32": LNPPreset(
        name="FO-32 (pulmonary, ML-designed)",
        ionizable_lipid="FO-32",
        helper_lipid="DOPE",
        cholesterol_pct=24.0,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.0,
        ionizable_mol_pct=60.0,
        helper_mol_pct=10.0,
        n_p_ratio=8.0,
        source="Witten et al. Nat Biotech 2025",
        notes="Top hit from ML-guided screening (>1.6M candidates); ferret-lung delivery.",
    ),
    "FO-35": LNPPreset(
        name="FO-35 (pulmonary, ML-designed)",
        ionizable_lipid="FO-35",
        helper_lipid="DOPE",
        cholesterol_pct=24.0,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.0,
        ionizable_mol_pct=60.0,
        helper_mol_pct=10.0,
        n_p_ratio=8.0,
        source="Witten et al. Nat Biotech 2025",
        notes="Second ML-discovered lipid; also ferret-lung active.",
    ),
    # saRNA-specific: usually lower ionizable%, higher cholesterol for stability
    "saRNA-gen-1": LNPPreset(
        name="saRNA generic",
        ionizable_lipid="SM-102",
        helper_lipid="DSPC",
        cholesterol_pct=46.0,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.0,
        ionizable_mol_pct=44.0,
        helper_mol_pct=9.0,
        n_p_ratio=6.0,
        source="Arcturus ARCT-154 disclosures / Pollock et al.",
        notes="Higher cholesterol + lower PEG for saRNA stability.",
    ),
    # Tumor-targeting via TME normalization (e.g. ionizable + losartan co-treatment)
    "tumor-it": LNPPreset(
        name="intratumoral TRAIL-mRNA LNP",
        ionizable_lipid="C12-200",
        helper_lipid="DOPE",
        cholesterol_pct=28.0,
        peg_lipid="DMG-PEG2000",
        peg_mol_pct=1.0,
        ionizable_mol_pct=55.0,
        helper_mol_pct=10.0,
        n_p_ratio=6.0,
        source="Costa et al. IJN 2025 (TRAIL mRNA)",
        notes="Pairs with losartan/Ang(1-7) TME normalization; intratumoral injection.",
    ),
}


# ---------- decision matrix -------------------------------------------------


@dataclass
class LNPAdvice:
    cargo: str
    target: str
    intent: str
    shortlist: list[dict]
    notes: list[str] = field(default_factory=list)


def recommend(
    *,
    target: str = "liver",
    cargo: str = "mRNA",
    intent: str = "cancer vaccine",
    n: int = 3,
) -> LNPAdvice:
    """Return a ranked shortlist of LNP formulations for a given scenario."""
    target = target.lower()
    cargo = cargo.lower()
    intent = intent.lower()
    notes: list[str] = []

    # ---- shortlist per target ----
    if target in {"lung", "pulmonary", "airway"}:
        pool = ["FO-32", "FO-35"]
        notes.append(
            "Pulmonary delivery benefits from ML-discovered biodegradable lipids (Witten 2025)."
        )
    elif target in {"liver", "hepatic"}:
        pool = ["C12-200", "SM-102", "ALC-0315"]
        notes.append("Liver is the natural sink for systemically administered LNPs.")
    elif target in {"spleen", "immune", "dendritic"}:
        pool = ["SM-102", "ALC-0315"]
        notes.append("Spleen/APCs are reached by ionizable LNPs at standard composition.")
    elif target in {"muscle", "im"}:
        pool = ["SM-102", "ALC-0315"]
        notes.append("Intramuscular LNPs typically use vaccine-grade composition.")
    elif target in {"tumor", "tme", "local"}:
        pool = ["tumor-it", "C12-200"]
        notes.append("Consider TME normalization (losartan/Ang(1-7)) alongside intratumoral LNP.")
    else:
        pool = ["SM-102", "ALC-0315", "C12-200"]
        notes.append(f"Unknown target '{target}'; defaulting to clinical-grade compositions.")

    # ---- cargo-specific adjustments ----
    if cargo in {"sarna", "srna", "self-amplifying"}:
        if "saRNA-gen-1" not in pool:
            pool.insert(0, "saRNA-gen-1")
        notes.append("saRNA prefers higher cholesterol and lower PEG for replicon stability.")
    elif cargo in {"circrna", "circular"}:
        notes.append("circRNA tolerates standard LNP composition; lower N/P often sufficient.")
    elif cargo in {"mrna", ""}:
        pass
    elif cargo in {"sirna", "sgrna"}:
        notes.append("siRNA/sgRNA cargoes typically use lower ionizable% (≤50%).")
    elif cargo in {"cas9", "crispr"}:
        if "C12-200" not in pool:
            pool.insert(0, "C12-200")
        notes.append("Cas9 mRNA + sgRNA co-delivery: C12-200 or similar hepatotropic LNPs.")
    else:
        notes.append(f"Unrecognized cargo '{cargo}'; using standard ratios.")

    # ---- intent-specific ----
    if intent in {"cancer vaccine", "vaccine"}:
        notes.append(
            "For cancer vaccines, MHC-I presentation of encoded antigen is the priority — clinical SM-102/ALC-0315 work well."
        )
    elif intent in {"gene editing", "crispr", "knock-in", "knock-out"}:
        if "C12-200" not in pool:
            pool.insert(0, "C12-200")
        notes.append(
            "Gene editing needs high transfection + minimal double-stranded RNA contaminants."
        )
    elif intent in {"protein replacement", "enzyme", "missing protein"}:
        notes.append("Protein-replacement mRNAs need high translation; codon-optimize + high dose.")

    shortlist = [PRESETS[k].as_dict() for k in pool if k in PRESETS][:n]
    return LNPAdvice(cargo=cargo, target=target, intent=intent, shortlist=shortlist, notes=notes)


# ---------- CLI ------------------------------------------------------------


def _run_cli(argv: list[str]) -> int:
    import argparse
    import json as _json
    from pathlib import Path as _Path

    p = argparse.ArgumentParser(prog="mrnavax lnp")
    p.add_argument(
        "--target", default="liver", help="liver | lung | spleen | muscle | tumor | dendritic"
    )
    p.add_argument("--cargo", default="mRNA", help="mRNA | saRNA | circRNA | siRNA | Cas9")
    p.add_argument("--intent", default="cancer vaccine")
    p.add_argument("--n", type=int, default=3)
    p.add_argument("--out")
    args = p.parse_args(argv)

    advice = recommend(target=args.target, cargo=args.cargo, intent=args.intent, n=args.n)
    out_text = _json.dumps(asdict(advice), indent=2)
    if args.out:
        _Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))


__all__ = [
    # Preset dataclass + curated LNP formulation table.
    "LNPPreset",
    "PRESETS",
    # Decision-matrix dataclass.
    "LNPAdvice",
    # Public API.
    "recommend",
]

