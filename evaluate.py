"""
evaluate.py
-----------
Tests the three load-bearing claims this system makes. Run it, screenshot the
output, and put that screenshot in the case study.

The brief asks you to test the load-bearing claims rather than the demo. These
are the claims:

  CLAIM 1  An absent required position produces a refusal, not a default call.
  CLAIM 2  A no-call ('--') is treated identically to an absent position.
  CLAIM 3  Strand representation does not change the call.

Claim 3 matters more than it looks. rs4149056 is reported T/C by CPIC and A/G
by some arrays. A system that assumed one strand would silently invert the
result for files using the other. This test proves the normalisation holds.

Run:  python evaluate.py
"""

from parser import parse_genotype_file, GenotypeFile
from caller import analyse

PASS = "PASS"
FAIL = "FAIL"
results = []


def check(name, condition, detail=""):
    status = PASS if condition else FAIL
    results.append((status, name, detail))
    print(f"[{status}] {name}")
    if detail:
        print(f"        {detail}")


def load(path):
    with open(path, "rb") as fh:
        return parse_genotype_file(fh.read(), path)


print("=" * 70)
print("EVALUATION: three load-bearing claims")
print("=" * 70)

# ---------------------------------------------------------------------------
print("\nCLAIM 1 — an absent required position produces a refusal")
print("-" * 70)

findings = {f.gene: f for f in analyse(load("sample_data/test_missing_star3.txt"))}
cyp = findings["CYP2C19"]

check(
    "CYP2C19 is declined when rs4986893 (*3) is absent",
    cyp.status == "declined",
    f"status={cyp.status}",
)
check(
    "The refusal names the specific missing position",
    "rs4986893" in (cyp.reason or ""),
    f"reason: {cyp.reason}",
)
check(
    "SLCO1B1 is still called — one gene's gap does not suppress the other",
    findings["SLCO1B1"].status == "called",
    f"SLCO1B1 = {findings['SLCO1B1'].diplotype} "
    f"({findings['SLCO1B1'].phenotype})",
)

# ---------------------------------------------------------------------------
print("\nCLAIM 2 — a no-call is treated identically to an absent position")
print("-" * 70)

findings = {f.gene: f for f in analyse(load("sample_data/test_nocall.txt"))}
slco = findings["SLCO1B1"]

check(
    "SLCO1B1 is declined when rs4149056 reads '--'",
    slco.status == "declined",
    f"status={slco.status}",
)
check(
    "A present-but-unreadable position is NOT defaulted to reference",
    slco.phenotype is None,
    "phenotype is None rather than 'Normal Function'",
)
check(
    "CYP2C19 is called normally in the same file",
    findings["CYP2C19"].status == "called",
    f"CYP2C19 = {findings['CYP2C19'].diplotype} "
    f"({findings['CYP2C19'].phenotype})",
)

# ---------------------------------------------------------------------------
print("\nCLAIM 3 — strand representation does not change the call")
print("-" * 70)

# Same person, two strand conventions. rs4149056 variant allele is C (or G on
# the opposite strand); rs12248560 variant allele is T (or A opposite).
forward = GenotypeFile(
    {"rs4244285": "AG", "rs4986893": "GG", "rs12248560": "CT", "rs4149056": "TC"},
    4, "forward strand",
)
reverse = GenotypeFile(
    {"rs4244285": "TC", "rs4986893": "CC", "rs12248560": "GA", "rs4149056": "AG"},
    4, "reverse strand",
)

fwd = {f.gene: f for f in analyse(forward)}
rev = {f.gene: f for f in analyse(reverse)}

check(
    "CYP2C19 phenotype identical across strand conventions",
    fwd["CYP2C19"].phenotype == rev["CYP2C19"].phenotype,
    f"forward={fwd['CYP2C19'].phenotype} | reverse={rev['CYP2C19'].phenotype}",
)
check(
    "SLCO1B1 phenotype identical across strand conventions",
    fwd["SLCO1B1"].phenotype == rev["SLCO1B1"].phenotype,
    f"forward={fwd['SLCO1B1'].phenotype} | reverse={rev['SLCO1B1'].phenotype}",
)

# ---------------------------------------------------------------------------
print("\nBONUS — contested combinations are flagged, not smoothed over")
print("-" * 70)

check(
    "One no-function plus one increased-function allele carries a caveat",
    fwd["CYP2C19"].contested is not None,
    f"diplotype={fwd['CYP2C19'].diplotype}, caveat present",
)

# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
failed = [r for r in results if r[0] == FAIL]
print(f"{len(results) - len(failed)}/{len(results)} checks passed")
if failed:
    print("\nFAILURES:")
    for _, name, detail in failed:
        print(f"  - {name}: {detail}")
print("=" * 70)
