"""
pgx_knowledge.py
----------------
The entire evidence base for this system, as a static table.

WHY THIS IS A FILE AND NOT A DATABASE OR AN API CALL:
Two gene-drug pairs is a small, stable, hand-curated dataset. Wrapping it in a
database or a live API would add moving parts without adding defensibility.
Every entry below carries the source it came from so a reviewer can check it.

WHY THERE IS NO LANGUAGE MODEL IN THIS FILE:
Mapping a genotype to a star allele to a metabolizer phenotype is a lookup with
published rules. A model would make it probabilistic for no gain and would
break traceability. The model-free design is a deliberate decision, not a gap.

VERIFY-BEFORE-PUBLISH:
Every `source_url` below should be opened and the claim confirmed by hand before
you ship. That spot-check IS the evaluation step for source discipline.
"""

# ---------------------------------------------------------------------------
# STRAND HANDLING
# ---------------------------------------------------------------------------
# Consumer genotyping arrays report alleles relative to a strand that is not
# always the one used in the pharmacogenomics literature. rs4149056 is the
# classic case: reported as T/C in CPIC documents, but as A/G by some arrays.
# We normalise by accepting both representations and mapping them to a single
# canonical allele. Silently guessing the strand is a known source of wrong
# calls, so this mapping is explicit and visible in the output.

STRAND_EQUIVALENTS = {
    "A": {"A", "T"},
    "T": {"A", "T"},
    "C": {"C", "G"},
    "G": {"C", "G"},
}


def alleles_match(observed: str, canonical: str) -> bool:
    """True if an observed base is the canonical base or its strand complement."""
    return observed.upper() in STRAND_EQUIVALENTS.get(canonical.upper(), set())


# ---------------------------------------------------------------------------
# VARIANT DEFINITIONS
# ---------------------------------------------------------------------------
# Each entry defines ONE variant that a star allele is built from.
#   canonical_variant : the base that indicates the variant (non-reference) allele
#   canonical_reference : the base that indicates the reference allele
#
# A gene's phenotype can only be called if EVERY rsid listed under
# `required_rsids` for that gene was actually genotyped in the uploaded file.

VARIANTS = {
    "rs4244285": {
        "gene": "CYP2C19",
        "star_allele": "*2",
        "canonical_reference": "G",
        "canonical_variant": "A",
        "function": "no function",
        "description": "c.681G>A, creates an aberrant splice site",
    },
    "rs4986893": {
        "gene": "CYP2C19",
        "star_allele": "*3",
        "canonical_reference": "G",
        "canonical_variant": "A",
        "function": "no function",
        "description": "c.636G>A, premature stop codon",
    },
    "rs12248560": {
        "gene": "CYP2C19",
        "star_allele": "*17",
        "canonical_reference": "C",
        "canonical_variant": "T",
        "function": "increased function",
        "description": "c.-806C>T, promoter variant raising transcription",
    },
    "rs4149056": {
        "gene": "SLCO1B1",
        "star_allele": "*5",
        "canonical_reference": "T",
        "canonical_variant": "C",
        "function": "decreased function",
        "description": "c.521T>C (Val174Ala), reduces hepatic transporter uptake",
    },
}


# ---------------------------------------------------------------------------
# GENE DEFINITIONS
# ---------------------------------------------------------------------------
# `required_rsids` is the coverage contract. If the uploaded file is missing any
# of these, the gene is DECLINED rather than defaulted to *1/*1.
#
# This is the central design decision of the project. Most consumer DNA
# interpretation tools treat an untested position as reference, which silently
# converts "we did not look" into "you are normal". This system refuses.
#
# `untested_alleles` names what this panel never interrogates at all, so even a
# successful call is reported with its own coverage ceiling attached.

GENES = {
    "CYP2C19": {
        "required_rsids": ["rs4244285", "rs4986893", "rs12248560"],
        "untested_alleles": [
            "*4", "*5", "*6", "*7", "*8", "*9", "*10", "*35",
        ],
        "coverage_note": (
            "This panel tests *2, *3 and *17 only. Rarer no-function alleles "
            "(*4 through *10, *35) are not interrogated, so a result reported "
            "here as Normal Metabolizer could still carry an untested "
            "no-function allele."
        ),
    },
    "SLCO1B1": {
        "required_rsids": ["rs4149056"],
        "untested_alleles": ["*15", "*17 haplotype context"],
        "coverage_note": (
            "rs4149056 is the single variant that drives the established "
            "simvastatin association. Haplotype context that distinguishes "
            "*5 from *15/*17 is not resolved by a single SNP."
        ),
    },
}


# ---------------------------------------------------------------------------
# PHENOTYPE RULES
# ---------------------------------------------------------------------------
# Keyed by (count of no-function alleles, count of increased-function alleles).
# Standard CPIC activity-score logic, simplified to the alleles this panel covers.

CYP2C19_PHENOTYPE = {
    (0, 0): "Normal Metabolizer",
    (0, 1): "Rapid Metabolizer",
    (0, 2): "Ultrarapid Metabolizer",
    (1, 0): "Intermediate Metabolizer",
    (1, 1): "Intermediate Metabolizer",  # phasing-dependent, see caller.py
    (2, 0): "Poor Metabolizer",
}

SLCO1B1_PHENOTYPE = {
    0: "Normal Function",
    1: "Decreased Function",
    2: "Poor Function",
}


# ---------------------------------------------------------------------------
# DRUG GUIDANCE
# ---------------------------------------------------------------------------
# `evidence_level`     : CPIC's own grading of the gene-drug pair
# `population_measured`: which ancestry groups the underlying evidence came from
# `allele_frequency`   : why the same finding is not equally informative for
#                        every user — the population dimension the brief asks for
# `guideline_positions`: separate entries per guideline body, so agreement and
#                        disagreement are both visible rather than merged

DRUG_GUIDANCE = {
    "CYP2C19": {
        "drug": "Clopidogrel (Plavix)",
        "mechanism": (
            "Clopidogrel is a prodrug. CYP2C19 performs the activation step, so "
            "reduced CYP2C19 function means less active metabolite and weaker "
            "antiplatelet effect."
        ),
        "evidence_level": "CPIC Level A",
        "evidence_level_meaning": (
            "Level A means CPIC judges the evidence sufficient to recommend a "
            "prescribing change. It is their highest tier."
        ),
        "source_name": "CPIC Guideline for Clopidogrel and CYP2C19",
        "source_url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
        "last_updated": "2022",
        "population_measured": (
            "Primary evidence from acute coronary syndrome / PCI cohorts, "
            "predominantly European and East Asian ancestry."
        ),
        "allele_frequency": {
            "*2": {
                "East Asian": "~29%",
                "European": "~15%",
                "African": "~17%",
                "note": "Common in every population studied.",
            },
            "*3": {
                "East Asian": "~5-9%",
                "European": "<1%",
                "African": "<1%",
                "note": (
                    "Near-absent outside East Asian ancestry. A panel that "
                    "skipped *3 would misclassify East Asian users far more "
                    "often than European ones."
                ),
            },
            "*17": {
                "European": "~21%",
                "African": "~16%",
                "East Asian": "~2%",
            },
        },
        "guideline_positions": [
            {
                "body": "CPIC",
                "position": (
                    "For Poor and Intermediate Metabolizers in ACS/PCI, "
                    "recommends an alternative antiplatelet such as prasugrel "
                    "or ticagrelor where not contraindicated."
                ),
                "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
            },
            {
                "body": "DPWG (Dutch Pharmacogenetics Working Group)",
                "position": (
                    "Also recommends an alternative for reduced-function "
                    "metabolizers, but scopes its recommendation by clinical "
                    "indication differently from CPIC. Confirm the current "
                    "wording before relying on this line."
                ),
                "url": "https://www.pharmgkb.org/page/dpwg",
            },
        ],
        "contested": (
            "Whether CYP2C19 genotype should guide clopidogrel therapy outside "
            "the ACS/PCI setting is not settled. Guideline bodies scope their "
            "recommendations differently, and this system reports that rather "
            "than picking a side."
        ),
    },
    "SLCO1B1": {
        "drug": "Simvastatin",
        "mechanism": (
            "SLCO1B1 encodes a transporter that moves statins into the liver. "
            "Reduced function leaves more drug in circulation, which raises "
            "muscle toxicity risk."
        ),
        "evidence_level": "CPIC Level A",
        "evidence_level_meaning": (
            "Level A means CPIC judges the evidence sufficient to recommend a "
            "prescribing change."
        ),
        "source_name": "CPIC Guideline for Statins and SLCO1B1, ABCG2, CYP2C9",
        "source_url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/",
        "last_updated": "2022",
        "population_measured": (
            "The landmark association came from the SEARCH trial, a "
            "predominantly European-ancestry cohort."
        ),
        "allele_frequency": {
            "*5 (rs4149056 C allele)": {
                "European": "~15%",
                "East Asian": "~12%",
                "African": "~2%",
                "note": (
                    "Markedly rarer in African ancestry. The absence of this "
                    "variant is correspondingly less reassuring for a user of "
                    "European ancestry than the raw result suggests."
                ),
            },
        },
        "guideline_positions": [
            {
                "body": "CPIC",
                "position": (
                    "For decreased or poor function, recommends limiting "
                    "simvastatin dose or choosing an alternative statin."
                ),
                "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/",
            },
            {
                "body": "DPWG",
                "position": (
                    "Issues dose-limit guidance for simvastatin in reduced "
                    "function carriers. Confirm the current wording before "
                    "relying on this line."
                ),
                "url": "https://www.pharmgkb.org/page/dpwg",
            },
        ],
        "contested": (
            "The strength of the association is dose-dependent and best "
            "established for simvastatin specifically. Extending it to other "
            "statins is weaker evidence."
        ),
    },
}


# ---------------------------------------------------------------------------
# MEDICAL SCOPE
# ---------------------------------------------------------------------------
DISCLAIMER = (
    "This is a sense-making tool, not a diagnostic. It reads a consumer "
    "genotyping file and shows what published pharmacogenomic guidelines say "
    "about a small number of variants, along with what it could not determine. "
    "It does not replace clinical testing or a clinician's judgment, and no "
    "medication decision should be made from it."
)
