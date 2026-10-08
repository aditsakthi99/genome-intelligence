"""
caller.py
---------
Turns a parsed genotype file into findings, declines, and the reasons for both.

THIS FILE IS THE PROJECT.
Everything else is plumbing. The rules below decide when the system is allowed
to tell the user something and when it must refuse, and every refusal names
what was missing.

THREE OUTCOMES PER GENE
  CALLED    - every required position was genotyped and the result resolves
  DECLINED  - a required position was absent or failed to read
  UNRESOLVED- positions were read but the result does not resolve to a diplotype

The distinction between DECLINED and UNRESOLVED matters. The first is a
coverage failure (the chip did not look). The second is an interpretation
failure (the chip looked, and the answer is ambiguous). Collapsing them would
hide which one happened.

WHY NO MODEL RUNS HERE
Star allele assignment is a published rule set with a finite input space. A
model would introduce variance into a deterministic step and break the
requirement that every claim trace to a source. Choosing not to use a model is
a decision, and this comment is where it is defended.
"""

from pgx_knowledge import (
    VARIANTS,
    GENES,
    CYP2C19_PHENOTYPE,
    SLCO1B1_PHENOTYPE,
    DRUG_GUIDANCE,
    alleles_match,
)


class Finding:
    """One gene's result, in whatever state it resolved to."""

    def __init__(self, gene, status, **kwargs):
        self.gene = gene
        self.status = status              # "called" | "declined" | "unresolved"
        self.diplotype = kwargs.get("diplotype")
        self.phenotype = kwargs.get("phenotype")
        self.observed = kwargs.get("observed", {})       # {rsid: genotype}
        self.missing = kwargs.get("missing", [])         # rsids not genotyped
        self.reason = kwargs.get("reason")               # why declined/unresolved
        self.contested = kwargs.get("contested")         # phenotype-level caveat
        self.guidance = DRUG_GUIDANCE.get(gene, {})
        self.coverage_note = GENES[gene]["coverage_note"]
        self.untested_alleles = GENES[gene]["untested_alleles"]


def _count_variant_alleles(genotype: str, rsid: str) -> int:
    """
    How many copies of the variant (non-reference) allele this genotype carries.

    Strand-aware: a C/G ambiguity is resolved by alleles_match rather than by
    assuming the vendor used the same strand as the literature.
    """
    definition = VARIANTS[rsid]
    variant_base = definition["canonical_variant"]
    return sum(1 for base in genotype if alleles_match(base, variant_base))


def _check_coverage(genotype_file, gene: str):
    """Return (observed, missing) for a gene's required positions."""
    observed, missing = {}, []
    for rsid in GENES[gene]["required_rsids"]:
        if genotype_file.has(rsid):
            observed[rsid] = genotype_file.get(rsid)
        else:
            missing.append(rsid)
    return observed, missing


def call_cyp2c19(genotype_file) -> Finding:
    """CYP2C19: three positions, combined into an activity-based phenotype."""
    observed, missing = _check_coverage(genotype_file, "CYP2C19")

    # ---- DECLINE: the chip did not look -----------------------------------
    if missing:
        missing_labels = [
            f"{rsid} ({VARIANTS[rsid]['star_allele']})" for rsid in missing
        ]
        return Finding(
            "CYP2C19",
            "declined",
            observed=observed,
            missing=missing,
            reason=(
                "Cannot call a CYP2C19 phenotype. This file does not contain a "
                "usable result for: " + ", ".join(missing_labels) + ". "
                "Reporting a phenotype anyway would require assuming these "
                "positions are reference, which this file does not show."
            ),
        )

    # ---- Count alleles ----------------------------------------------------
    no_function = (
        _count_variant_alleles(observed["rs4244285"], "rs4244285")   # *2
        + _count_variant_alleles(observed["rs4986893"], "rs4986893")  # *3
    )
    increased = _count_variant_alleles(observed["rs12248560"], "rs12248560")  # *17

    # ---- UNRESOLVED: looked, but the answer does not add up ---------------
    if no_function + increased > 2:
        return Finding(
            "CYP2C19",
            "unresolved",
            observed=observed,
            reason=(
                f"Observed {no_function} no-function and {increased} "
                "increased-function alleles, which is more than the two a "
                "person carries. This usually means a rare haplotype where "
                "two variants sit on the same chromosome. Resolving it "
                "requires phased data this file does not provide."
            ),
        )

    phenotype = CYP2C19_PHENOTYPE.get((no_function, increased))
    if phenotype is None:
        return Finding(
            "CYP2C19",
            "unresolved",
            observed=observed,
            reason="Allele combination does not map to a defined phenotype.",
        )

    # ---- Build a readable diplotype --------------------------------------
    alleles = []
    for rsid, star in (("rs4244285", "*2"), ("rs4986893", "*3"), ("rs12248560", "*17")):
        alleles.extend([star] * _count_variant_alleles(observed[rsid], rsid))
    while len(alleles) < 2:
        alleles.insert(0, "*1")
    # Sort numerically, not alphabetically: *2 comes before *17.
    # Reference (*1) is conventionally written first.
    diplotype = "/".join(
        sorted(alleles, key=lambda a: (a != "*1", int(a.lstrip("*"))))
    )

    # ---- Flag the genuinely contested combination -------------------------
    contested = None
    if no_function == 1 and increased == 1:
        contested = (
            "This diplotype carries one no-function and one increased-function "
            "allele. How the two offset each other is not settled, and the "
            "Intermediate Metabolizer assignment is the conservative reading "
            "rather than a consensus one."
        )

    return Finding(
        "CYP2C19",
        "called",
        diplotype=diplotype,
        phenotype=phenotype,
        observed=observed,
        contested=contested,
    )


def call_slco1b1(genotype_file) -> Finding:
    """SLCO1B1: one position, directly mapped to a transporter function class."""
    observed, missing = _check_coverage(genotype_file, "SLCO1B1")

    if missing:
        return Finding(
            "SLCO1B1",
            "declined",
            observed=observed,
            missing=missing,
            reason=(
                "Cannot call SLCO1B1 function. This file does not contain a "
                "usable result for rs4149056 (*5), which is the single variant "
                "the simvastatin association rests on."
            ),
        )

    variant_count = _count_variant_alleles(observed["rs4149056"], "rs4149056")
    phenotype = SLCO1B1_PHENOTYPE[variant_count]
    diplotype = {0: "*1/*1", 1: "*1/*5", 2: "*5/*5"}[variant_count]

    return Finding(
        "SLCO1B1",
        "called",
        diplotype=diplotype,
        phenotype=phenotype,
        observed=observed,
    )


def analyse(genotype_file) -> list[Finding]:
    """Run every gene in the panel. Order is stable so output is reproducible."""
    return [call_cyp2c19(genotype_file), call_slco1b1(genotype_file)]
