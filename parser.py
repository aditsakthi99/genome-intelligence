"""
parser.py
---------
Reads a raw consumer genotype file and returns a dict of {rsid: genotype}.

WHY THIS IS DELIBERATELY DUMB:
The parser's only job is to report what the file literally contains. It does
not fill gaps, infer missing positions, or normalise anything beyond stripping
whitespace. Every judgment call happens downstream in caller.py, where it can
be explained. A parser that quietly "helps" is how untested positions turn into
confident results.

SUPPORTED FORMATS
23andMe, AncestryDNA, and PGP exports are all tab-separated text with comment
lines starting with '#'. Column layouts differ slightly, so we detect them.

    23andMe / PGP : rsid  chromosome  position  genotype      (4 columns)
    AncestryDNA   : rsid  chromosome  position  allele1  allele2  (5 columns)

NO-CALLS
A position the chip attempted but failed to read is written as '--', '00', or
'NN' depending on vendor. These are treated as NOT GENOTYPED, identical to the
position being absent from the file. That equivalence is the point: both mean
"we do not know", and both must trigger a decline.
"""

NO_CALL_TOKENS = {"--", "00", "NN", "", "-", "I", "D", "DD", "II", "DI", "ID"}


class GenotypeFile:
    """Holds the parsed contents of one raw genotype file."""

    def __init__(self, genotypes: dict, total_lines: int, source_name: str):
        self.genotypes = genotypes          # {rsid: "AG"} — called positions only
        self.total_lines = total_lines      # how many data rows the file had
        self.source_name = source_name      # filename, for provenance display

    def has(self, rsid: str) -> bool:
        """True only if this rsid was present AND successfully called."""
        return rsid in self.genotypes

    def get(self, rsid: str):
        return self.genotypes.get(rsid)


def _clean_genotype(raw: str) -> str | None:
    """Return a normalised 2-base genotype, or None if this is a no-call."""
    g = raw.strip().upper().replace(" ", "")
    if g in NO_CALL_TOKENS:
        return None
    # Keep only A/C/G/T characters; anything else (indels, CNV markers) is not
    # something this system can reason about, so it is treated as uncalled.
    if not all(base in "ACGT" for base in g):
        return None
    if len(g) not in (1, 2):
        return None
    return g


def parse_genotype_file(file_bytes: bytes, source_name: str = "uploaded file") -> GenotypeFile:
    """
    Parse raw bytes from an uploaded genotype file.

    Returns a GenotypeFile. Raises ValueError if the file has no parseable rows,
    which usually means the user uploaded the wrong thing entirely.
    """
    try:
        text = file_bytes.decode("utf-8", errors="replace")
    except Exception as exc:
        raise ValueError(f"Could not read file as text: {exc}")

    genotypes = {}
    data_rows = 0

    for line in text.splitlines():
        line = line.strip()

        # Skip comments, blank lines, and the header row some vendors include.
        if not line or line.startswith("#"):
            continue
        if line.lower().startswith("rsid"):
            continue

        # Vendors use tabs, but a few exports use commas.
        parts = line.split("\t") if "\t" in line else line.split(",")
        parts = [p.strip().strip('"') for p in parts]

        if len(parts) < 4:
            continue

        rsid = parts[0]
        if not rsid.lower().startswith("rs"):
            # Internal vendor identifiers (i-prefixed on 23andMe) are skipped;
            # nothing in our knowledge base is keyed to them.
            continue

        data_rows += 1

        if len(parts) >= 5:
            raw_genotype = parts[3] + parts[4]   # AncestryDNA layout
        else:
            raw_genotype = parts[3]              # 23andMe / PGP layout

        cleaned = _clean_genotype(raw_genotype)
        if cleaned is not None:
            genotypes[rsid] = cleaned

    if data_rows == 0:
        raise ValueError(
            "No genotype rows found. Expected a raw DNA export from 23andMe, "
            "AncestryDNA, MyHeritage, or the Personal Genome Project."
        )

    return GenotypeFile(genotypes, data_rows, source_name)
