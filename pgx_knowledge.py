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
        "source_name": (
            "CPIC Guideline for CYP2C19 Genotype and Clopidogrel Therapy: "
            "2022 Update"
        ),
        "source_citation": (
            "Lee CR, Luzum JA, Sangkuhl K, et al. Clin Pharmacol Ther. "
            "2022;112(5):959-967. doi:10.1002/cpt.2526"
        ),
        "source_url": "https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf",
        "last_updated": "2022",
        "verified_on": "2026-10-08, against the CPIC 2022 publication PDF",
        "population_measured": (
            "Primary evidence from cardiovascular (ACS/PCI) and neurovascular "
            "cohorts, predominantly European and East Asian ancestry. CPIC "
            "also states plainly that CYP2C19 variation does not account for "
            "all variability in clopidogrel response, so genotype is one input "
            "among several."
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
                "body": "CPIC (cardiovascular indications)",
                "position": (
                    "For Poor and Intermediate Metabolizers, avoid standard "
                    "dose clopidogrel if possible; use prasugrel or ticagrelor "
                    "where not contraindicated. The 2022 update raised the "
                    "strength of the Intermediate Metabolizer recommendation "
                    "to strong."
                ),
                "url": "https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf",
            },
            {
                "body": "CPIC (neurovascular indications)",
                "position": (
                    "Separate recommendation table: an alternative P2Y12 "
                    "inhibitor where clinically indicated and not "
                    "contraindicated. The 2022 update expanded CPIC's scope "
                    "beyond the cardiovascular setting."
                ),
                "url": "https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf",
            },
            {
                "body": "DPWG (Dutch Pharmacogenetics Working Group)",
                "position": (
                    "Also publishes a CYP2C19-clopidogrel guideline and also "
                    "recommends an alternative drug for Intermediate and Poor "
                    "Metabolizers. The two bodies agree on direction here; "
                    "they differ in how recommendations are scoped and graded."
                ),
                "url": "https://www.clinpgx.org/guidelineAnnotation/PA166104948/annotation",
            },
        ],
        "contested": (
            "The *2/*17 combination is the weak point. CPIC classifies it as "
            "an Intermediate Metabolizer but has carried that assignment as a "
            "provisional classification, because how one no-function and one "
            "increased-function allele offset each other is not established. "
            "This system reports the provisional status rather than hiding it "
            "behind a clean label."
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
        "source_name": (
            "CPIC Guideline for SLCO1B1, ABCG2, and CYP2C9 Genotypes and "
            "Statin-Associated Musculoskeletal Symptoms"
        ),
        "source_citation": (
            "Cooper-DeHoff RM, Niemi M, Ramsey LB, et al. Clin Pharmacol Ther. "
            "2022;111(5):1007-1021. doi:10.1002/cpt.2557"
        ),
        "source_url": "https://files.cpicpgx.org/data/guideline/publication/statins/2022/publication.pdf",
        "last_updated": "2022 (replaces the 2012 and 2014 simvastatin guidelines)",
        "verified_on": "2026-10-08, against the CPIC 2022 publication and supplement",
        "population_measured": (
            "The landmark association came from the SEARCH trial, in which "
            "myopathy was observed on 80 mg simvastatin in a predominantly "
            "European-ancestry cohort. CPIC's own 2022 supplement notes that "
            "rs4149056 does not appear to contribute to statin sensitivity in "
            "Asian populations in the same way, so this finding does not "
            "transfer cleanly across ancestries."
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
                "body": "DPWG (Dutch Pharmacogenetics Working Group)",
                "position": (
                    "Also publishes an SLCO1B1-simvastatin guideline, "
                    "recommending a lower dose, an alternative drug, or "
                    "additional monitoring such as routine CK surveillance. "
                    "Broadly aligned with CPIC in direction."
                ),
                "url": "https://www.clinpgx.org/guidelineAnnotation/PA166105005/annotation",
            },
        ],
        "contested": (
            "The association is dose-dependent and statin-specific. CPIC rates "
            "the simvastatin evidence high, atorvastatin moderate, and finds "
            "little evidence for pravastatin or rosuvastatin, so this is not a "
            "class effect. Risk also scales sharply with dose: roughly 2.6-fold "
            "per C allele at 40 mg versus substantially higher at 80 mg. A "
            "result reported without the dose attached means less than it "
            "appears to."
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
