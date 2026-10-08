"""
make_test_files.py
------------------
Generates three synthetic genotype files in 23andMe export format.

WHY SYNTHETIC FILES EXIST ALONGSIDE REAL ONES:
Real PGP files are the honest end-to-end test and you should use one for your
user session. But you cannot choose what a real file contains, so you cannot
use it to prove the decline path fires. These three files are constructed so
that each one exercises exactly one behaviour, which is what makes the
evaluation checkable rather than anecdotal.

    test_full_coverage.txt  -> all four positions called      -> expect 2 findings
    test_missing_star3.txt  -> rs4986893 absent from the file -> expect CYP2C19 DECLINED
    test_nocall.txt         -> rs4149056 present but '--'     -> expect SLCO1B1 DECLINED

The third file is the important one. A position that is present-but-unreadable
and a position that is absent must produce the same refusal, because both mean
the same thing: nobody looked successfully.

Run:  python make_test_files.py
"""

import os
import random

# Create the output folder if it does not exist, so this runs on a fresh clone.
os.makedirs("sample_data", exist_ok=True)

random.seed(42)  # reproducible filler, so test files are stable across runs

HEADER = """# This data file generated for testing purposes.
# Format mirrors a 23andMe raw data export.
# rsid\tchromosome\tposition\tgenotype
"""

# The positions this system actually reads.
TARGETS = {
    "rs4244285":  ("10", "94781859"),   # CYP2C19 *2
    "rs4986893":  ("10", "94780653"),   # CYP2C19 *3
    "rs12248560": ("10", "94761900"),   # CYP2C19 *17
    "rs4149056":  ("12", "21178615"),   # SLCO1B1 *5
}


def filler_rows(n=300):
    """Irrelevant rows, so the parser is exercised on a realistically sized file."""
    rows = []
    for i in range(n):
        rsid = f"rs{random.randint(1000000, 9999999)}"
        chrom = str(random.randint(1, 22))
        pos = str(random.randint(100000, 99999999))
        geno = random.choice(["AA", "AG", "GG", "CC", "CT", "TT", "--"])
        rows.append(f"{rsid}\t{chrom}\t{pos}\t{geno}")
    return rows


def write_file(path, target_genotypes, omit=()):
    """
    target_genotypes : {rsid: genotype} for the positions we care about
    omit             : rsids to leave out of the file entirely
    """
    rows = filler_rows(150)
    for rsid, (chrom, pos) in TARGETS.items():
        if rsid in omit:
            continue
        rows.append(f"{rsid}\t{chrom}\t{pos}\t{target_genotypes[rsid]}")
    rows.extend(filler_rows(150))

    with open(path, "w") as fh:
        fh.write(HEADER)
        fh.write("\n".join(rows) + "\n")
    print(f"wrote {path}")


# --- 1. Full coverage: heterozygous *2, reference *3, heterozygous *17 ------
# Expected: CYP2C19 called as *2/*17 Intermediate Metabolizer WITH the contested
# flag set, and SLCO1B1 called as *1/*5 Decreased Function.
write_file(
    "sample_data/test_full_coverage.txt",
    {
        "rs4244285": "AG",    # one *2 copy
        "rs4986893": "GG",    # no *3
        "rs12248560": "CT",   # one *17 copy
        "rs4149056": "TC",    # one *5 copy
    },
)

# --- 2. Missing *3: the position is simply not in the file ------------------
# Expected: CYP2C19 DECLINED naming rs4986893. SLCO1B1 still called.
# This is the narrative centrepiece: *3 is near-absent in European ancestry but
# reaches 5-9% in East Asian ancestry, so defaulting it to reference would
# misclassify East Asian users far more often.
write_file(
    "sample_data/test_missing_star3.txt",
    {
        "rs4244285": "AA",    # two *2 copies
        "rs4986893": "GG",    # omitted below
        "rs12248560": "CC",
        "rs4149056": "TT",
    },
    omit=("rs4986893",),
)

# --- 3. No-call: the position is present but the chip failed to read it -----
# Expected: SLCO1B1 DECLINED. CYP2C19 called as *1/*1 Normal Metabolizer.
write_file(
    "sample_data/test_nocall.txt",
    {
        "rs4244285": "GG",
        "rs4986893": "GG",
        "rs12248560": "CC",
        "rs4149056": "--",    # chip attempted, failed
    },
)

print("\nDone. Three files in sample_data/")
