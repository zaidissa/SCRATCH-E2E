#!/usr/bin/env python3
"""
Regression test for bin/cellranger_renaming.py.

    python3 tests/test_cellranger_renaming.py            # tests ../bin/cellranger_renaming.py
    python3 tests/test_cellranger_renaming.py <script>   # tests another copy

No dependencies beyond the standard library. Each case writes empty-ish FASTQs
whose CONTENT is their original name, runs the script, and checks that every
input survives under the name cellranger expects -- or that the script stops
with a non-zero exit and leaves the inputs unharmed.

Why it exists: earlier versions renamed with os.rename() onto names that were
already taken, so one read file silently replaced another and the script still
exited 0. Seen in two forms:
  * read number taken from the last character before ".fastq.gz": Illumina
    names ("..._R2_001.fastq.gz") end in "1", so R2 became R1, and every lane
    was forced to L001;
  * read number taken only from "_R1_"/"_R2_" with a default of 1: "X_R2.fastq.gz"
    and fetchngs "SRR_2.fastq.gz" became R1.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = (sys.argv[1] if len(sys.argv) > 1
          else os.path.join(HERE, "..", "bin", "cellranger_renaming.py"))

# (description, sample, input files, expected {output name: original name} or None = must exit non-zero)
CASES = [
    ("already 10x, prefix = sample", "S",
     ["S_S1_L001_R1_001.fastq.gz", "S_S1_L001_R2_001.fastq.gz"],
     {"S_S1_L001_R1_001.fastq.gz": "S_S1_L001_R1_001.fastq.gz",
      "S_S1_L001_R2_001.fastq.gz": "S_S1_L001_R2_001.fastq.gz"}),
    ("already 10x, four lanes", "S",
     [f"S_S6_L00{l}_R{r}_001.fastq.gz" for l in (1, 2, 3, 4) for r in (1, 2)],
     {f"S_S6_L00{l}_R{r}_001.fastq.gz": f"S_S6_L00{l}_R{r}_001.fastq.gz"
      for l in (1, 2, 3, 4) for r in (1, 2)}),
    ("Illumina names, prefix != sample", "S",
     ["lib_S3_L002_R1_001.fastq.gz", "lib_S3_L002_R2_001.fastq.gz"],
     {"S_S3_L002_R1_001.fastq.gz": "lib_S3_L002_R1_001.fastq.gz",
      "S_S3_L002_R2_001.fastq.gz": "lib_S3_L002_R2_001.fastq.gz"}),
    ("Illumina names, prefix != sample, two lanes", "S",
     [f"a_S1_L00{l}_R{r}_001.fastq.gz" for l in (1, 2) for r in (1, 2)],
     {f"S_S1_L00{l}_R{r}_001.fastq.gz": f"a_S1_L00{l}_R{r}_001.fastq.gz"
      for l in (1, 2) for r in (1, 2)}),
    ("X_R1 / X_R2, no lane", "S",
     ["X_R1.fastq.gz", "X_R2.fastq.gz"],
     {"S_S1_L001_R1_001.fastq.gz": "X_R1.fastq.gz",
      "S_S1_L001_R2_001.fastq.gz": "X_R2.fastq.gz"}),
    ("fetchngs SRR_1 / SRR_2", "S",
     ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"],
     {"S_S1_L001_R1_001.fastq.gz": "SRR1_1.fastq.gz",
      "S_S1_L001_R2_001.fastq.gz": "SRR1_2.fastq.gz"}),
    ("sample name itself contains _R1_", "S",
     ["pt_R1_tcr_S1_L001_R1_001.fastq.gz", "pt_R1_tcr_S1_L001_R2_001.fastq.gz"],
     {"S_S1_L001_R1_001.fastq.gz": "pt_R1_tcr_S1_L001_R1_001.fastq.gz",
      "S_S1_L001_R2_001.fastq.gz": "pt_R1_tcr_S1_L001_R2_001.fastq.gz"}),
    ("two runs, no lane numbers -> must refuse", "S",
     ["A_R1.fastq.gz", "A_R2.fastq.gz", "B_R1.fastq.gz", "B_R2.fastq.gz"], None),
    ("no read number -> must refuse", "S",
     ["mystery.fastq.gz"], None),
]


def run_case(sample, files):
    d = tempfile.mkdtemp()
    for f in files:
        with open(os.path.join(d, f), "w") as fh:
            fh.write(f)
    p = subprocess.run([sys.executable, "-I", SCRIPT, sample, d],
                       capture_output=True, text=True)
    out = {}
    for f in os.listdir(d):
        with open(os.path.join(d, f)) as fh:
            out[f] = fh.read()
    return p.returncode, out


def main():
    failures = 0
    for desc, sample, files, expected in CASES:
        code, out = run_case(sample, files)
        survivors = sorted(out.values())
        if expected is None:
            ok = code != 0 and survivors == sorted(files)
            want = "non-zero exit, every input still present"
        else:
            ok = code == 0 and out == expected
            want = "exit 0, every input under its 10x name"
        lost = len(files) - len(out)
        print(f"{'PASS' if ok else 'FAIL'}  {desc}"
              + ("" if ok else f"  (wanted {want}; got exit {code}, "
                               f"{len(out)}/{len(files)} files{', %d LOST' % lost if lost else ''})"))
        failures += not ok
    print(f"\n{len(CASES) - failures}/{len(CASES)} passed  [{SCRIPT}]")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
