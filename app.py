"""
app.py
------
Streamlit interface.

ONE LAYOUT DECISION WORTH DEFENDING:
Declines are rendered at the same visual weight as findings, in the same
column, above the fold. They are not collapsed into an appendix. The entire
argument of this system is that "we could not determine this" is a result, and
burying it in an expander would contradict that in the interface while claiming
it in the copy.

Run locally:  streamlit run app.py
"""

import streamlit as st

from parser import parse_genotype_file
from caller import analyse
from pgx_knowledge import VARIANTS, DISCLAIMER

st.set_page_config(
    page_title="Genome Intelligence — Pharmacogenomics",
    page_icon="🧬",
    layout="centered",
)

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("🧬 Genome Intelligence")
st.markdown(
    "**A pharmacogenomics reader that shows you what it could not determine.**"
)
st.markdown(
    "Most consumer DNA interpretation tools treat a position the chip never "
    "tested as if it were normal. That silently converts *we did not look* "
    "into *you are fine*. This one refuses to make that call and tells you "
    "exactly what was missing."
)

st.info(DISCLAIMER, icon="⚕️")

# ---------------------------------------------------------------------------
# Upload
# ---------------------------------------------------------------------------
st.markdown("---")
st.subheader("1. Upload a raw genotype file")
st.caption(
    "Accepts raw DNA exports from 23andMe, AncestryDNA, MyHeritage, or the "
    "Personal Genome Project. The file is parsed in memory and is not stored."
)

uploaded = st.file_uploader("Raw genotype file (.txt or .csv)", type=["txt", "csv"])

with st.expander("No file handy? Use a test file"):
    st.markdown(
        """
        Three test files ship with this repository in `sample_data/`:

        | File | What it demonstrates |
        |---|---|
        | `test_full_coverage.txt` | Both genes called, one with a contested diplotype |
        | `test_missing_star3.txt` | CYP2C19 declined — a required position is absent |
        | `test_nocall.txt` | SLCO1B1 declined — position present but unreadable |

        For real data, the Harvard Personal Genome Project publishes consented
        genotype files at `my.pgp-hms.org/public_genetic_data`.
        """
    )

if uploaded is None:
    st.stop()

# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------
try:
    genotype_file = parse_genotype_file(uploaded.getvalue(), uploaded.name)
except ValueError as exc:
    st.error(f"Could not read that file: {exc}")
    st.stop()

st.success(
    f"Parsed **{uploaded.name}** — {genotype_file.total_lines:,} positions in "
    f"the file, {len(genotype_file.genotypes):,} with a usable call."
)

findings = analyse(genotype_file)
called = [f for f in findings if f.status == "called"]
refused = [f for f in findings if f.status in ("declined", "unresolved")]

st.markdown("---")
st.subheader("2. Results")

col_a, col_b = st.columns(2)
col_a.metric("Findings reported", len(called))
col_b.metric("Findings declined", len(refused))


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------
def render_observed(finding):
    """Show the raw genotypes this call was built from — the audit trail."""
    if not finding.observed:
        st.caption("No relevant positions were readable in this file.")
        return
    rows = []
    for rsid, genotype in finding.observed.items():
        definition = VARIANTS[rsid]
        rows.append(
            f"| `{rsid}` | {definition['star_allele']} | `{genotype}` | "
            f"{definition['description']} |"
        )
    st.markdown(
        "| Position | Allele | Your genotype | What it is |\n"
        "|---|---|---|---|\n" + "\n".join(rows)
    )


def render_evidence(finding):
    """Strength, provenance, recency and population — the four the brief asks for."""
    g = finding.guidance
    if not g:
        return

    st.markdown(f"**Drug:** {g['drug']}")
    st.markdown(f"*{g['mechanism']}*")

    st.markdown("**Evidence strength**")
    st.markdown(f"- {g['evidence_level']} — {g['evidence_level_meaning']}")
    st.markdown(
        f"- Source: [{g['source_name']}]({g['source_url']}) "
        f"(last updated {g['last_updated']})"
    )
    if g.get("source_citation"):
        st.caption(f"  {g['source_citation']}")
    if g.get("verified_on"):
        st.caption(f"  Claims in this entry hand-checked {g['verified_on']}.")

    st.markdown("**Population the evidence was measured in**")
    st.markdown(f"- {g['population_measured']}")
    for allele, freqs in g["allele_frequency"].items():
        parts = [f"{k} {v}" for k, v in freqs.items() if k != "note"]
        st.markdown(f"- `{allele}` frequency: {', '.join(parts)}")
        if "note" in freqs:
            st.caption(f"  {freqs['note']}")

    st.markdown("**What the guideline bodies say**")
    for position in g["guideline_positions"]:
        st.markdown(
            f"- **{position['body']}**: {position['position']} "
            f"([source]({position['url']}))"
        )

    st.warning(f"**Where this is contested:** {g['contested']}", icon="⚖️")


def render_coverage_ceiling(finding):
    st.error(
        f"**What this result cannot rule out:** {finding.coverage_note}",
        icon="🚧",
    )


# ---------------------------------------------------------------------------
# Declines first — the thesis of the system
# ---------------------------------------------------------------------------
if refused:
    st.markdown("### ⛔ Declined")
    st.caption(
        "These are not errors. The system had insufficient data and is saying "
        "so rather than defaulting to a normal result."
    )
    for finding in refused:
        with st.container(border=True):
            label = "Could not determine" if finding.status == "declined" \
                else "Could not resolve"
            st.markdown(f"#### {finding.gene} — {label}")
            st.markdown(finding.reason)

            if finding.missing:
                st.markdown("**Positions not usable in this file:**")
                for rsid in finding.missing:
                    d = VARIANTS[rsid]
                    st.markdown(
                        f"- `{rsid}` — {d['star_allele']} ({d['function']}), "
                        f"{d['description']}"
                    )

            if finding.observed:
                with st.expander("Positions that were readable"):
                    render_observed(finding)

            with st.expander("What this gene would have told you"):
                render_evidence(finding)

# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------
if called:
    st.markdown("### ✅ Reported")
    for finding in called:
        with st.container(border=True):
            st.markdown(f"#### {finding.gene} — {finding.phenotype}")
            st.markdown(f"**Diplotype:** `{finding.diplotype}`")

            if finding.contested:
                st.warning(
                    f"**This call is not clean:** {finding.contested}", icon="⚠️"
                )

            st.markdown("**What this was built from**")
            render_observed(finding)

            render_evidence(finding)
            render_coverage_ceiling(finding)

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.markdown("---")
st.caption(
    "Built for the 100xEngineers Cohort 7 capstone. No language model runs in "
    "the interpretation path: star allele assignment is a published rule set, "
    "and making it probabilistic would break the requirement that every claim "
    "trace back to a named source."
)
