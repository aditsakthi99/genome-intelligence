# Genome Intelligence — A Pharmacogenomics Reader That Declines

**100xEngineers Cohort 7 Capstone · Genome Intelligence brief · Solo build**

Most consumer DNA interpretation tools treat a position the genotyping chip
never tested as if it were the reference allele. That silently converts *we did
not look* into *you are normal*. This system refuses to make that call, and
names what was missing instead.

---

## The problem

A consumer genotyping array tests a subset of the variants that define a
pharmacogene. CYP2C19 has dozens of known star alleles; a typical consumer chip
interrogates three. The standard pipeline assigns `*1` — normal function — to
every allele it did not test, and reports a confident phenotype built partly on
positions nobody looked at.

The cost is not evenly distributed. CYP2C19 `*3` is near-absent in European
ancestry but reaches 5–9% in East Asian ancestry. A panel that omits `*3` and
defaults it to reference will misclassify East Asian users far more often than
European ones, and nothing in the output signals that.

## What this does instead

Three outcomes per gene, kept distinct:

| Outcome | Meaning |
|---|---|
| **Called** | Every required position was genotyped and the result resolves |
| **Declined** | A required position was absent, or the chip failed to read it |
| **Unresolved** | Positions were read, but the result does not resolve to a diplotype |

`Declined` and `Unresolved` are separated on purpose. The first is a coverage
failure — the chip did not look. The second is an interpretation failure — the
chip looked and the answer is ambiguous. Collapsing them would hide which one
happened.

Every reported finding carries four things the brief asks for: **evidence
strength** (CPIC level), **provenance** (named source, linked), **recency**
(guideline date), and **population** (which ancestry groups the underlying
evidence was measured in, plus allele frequencies across groups).

Each finding also carries its own coverage ceiling — what a "normal" result
still cannot rule out, given what this panel does not test.

## Scope

Two gene-drug pairs, deliberately.

- **CYP2C19 → clopidogrel** — `*2` (rs4244285), `*3` (rs4986893), `*17` (rs12248560)
- **SLCO1B1 → simvastatin** — `*5` (rs4149056)

CYP2D6 is excluded. Copy-number variation and hybrid alleles make it largely
invisible to array genotyping, and a system built on refusing to overclaim
should not start by overclaiming the hardest gene in the panel.

## Why no model runs in the interpretation path

Star allele assignment is a published rule set with a finite input space. A
language model would introduce variance into a deterministic step and break the
requirement that every claim trace to a named source. Choosing not to use one
is a design decision, defended in `caller.py`.

## Architecture

```
parser.py         Raw genotype file -> {rsid: genotype}. Deliberately dumb.
                  Reports what the file literally contains, fills no gaps.

pgx_knowledge.py  The entire evidence base as a static table. Every entry
                  carries the source URL it came from.

caller.py         The project. Coverage checks, diplotype assignment,
                  phenotype rules, and the decline logic.

app.py            Streamlit interface. Declines render at the same visual
                  weight as findings, above the fold.

evaluate.py       Tests the three load-bearing claims.
```

### Strand handling

rs4149056 is reported `T/C` in CPIC documents and `A/G` by some arrays. A system
that assumed one strand would silently invert the result for files using the
other. Normalisation is explicit in `pgx_knowledge.STRAND_EQUIVALENTS` and
tested in `evaluate.py`.

## Evaluation

`python evaluate.py` tests three claims:

1. An absent required position produces a refusal, not a default call.
2. A no-call (`--`) is treated identically to an absent position.
3. Strand representation does not change the call.

All checks pass. Output is in the case study.

## Running it

```bash
pip install -r requirements.txt
python make_test_files.py     # generates sample_data/
python evaluate.py            # runs the test suite
streamlit run app.py
```

## Test data

Three synthetic files in `sample_data/`, each exercising one behaviour:

| File | Demonstrates |
|---|---|
| `test_full_coverage.txt` | Both genes called, one contested diplotype |
| `test_missing_star3.txt` | CYP2C19 declined — required position absent |
| `test_nocall.txt` | SLCO1B1 declined — position present but unreadable |

Real data: the Harvard Personal Genome Project publishes consented, CC0-licensed
genotype files at `my.pgp-hms.org/public_genetic_data`.

## Scope and limits

This is a sense-making tool, not a diagnostic. It does not replace clinical
pharmacogenomic testing or a clinician's judgment, and no medication decision
should be made from it. Consumer array data is not validated for clinical use.

## Sources

- CPIC Guideline for Clopidogrel and CYP2C19 — https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/
- CPIC Guideline for Statins — https://cpicpgx.org/guidelines/cpic-guideline-for-statins/
- PharmGKB / DPWG — https://www.pharmgkb.org/page/dpwg
- Harvard Personal Genome Project — https://pgp.med.harvard.edu/
