#!/usr/bin/env python3
"""
auto_demux.py

Auto-detects sample barcodes from I1 index FASTQs and demultiplexes R1/R2
reads into per-sample FASTQ files.  No index samplesheet is required.

Multi-lane support
------------------
Pass one file per lane to --r1, --r2, --i1 (space-separated).  All lanes
must be listed in the same order so that R1[n] / R2[n] / I1[n] correspond
to the same lane.  Detection pools I1 reads across all lanes; demultiplexing
streams every lane sequentially and writes into the same per-sample output
files, producing one merged FASTQ per sample regardless of lane count.

Detection logic
---------------
1. Sample up to --n-detect reads total from I1 (spread evenly across lanes).
2. Build a frequency table of all I1 barcode sequences.
3. Any barcode above --min-freq (default 1%) that is not within 1 mismatch
   of an already-chosen centroid is called as a distinct sample.
4. Samples are auto-named <run_id>_S01_<centroid>, S02_... in descending
   frequency order.

Output FASTQs follow the cellranger naming convention:
  <sample>_S1_L001_R1_001.fastq.gz / <sample>_S1_L001_R2_001.fastq.gz
so they feed directly into cellranger count/vdj with --fastqs=. --sample=.

Usage (single lane)
-------------------
    auto_demux.py --run-id RUN1 --r1 R1.fastq.gz --r2 R2.fastq.gz
                  --i1 I1.fastq.gz --modality GEX

Usage (multiple lanes)
----------------------
    auto_demux.py --run-id RUN1
                  --r1 L001_R1.fastq.gz L002_R1.fastq.gz ...
                  --r2 L001_R2.fastq.gz L002_R2.fastq.gz ...
                  --i1 L001_I1.fastq.gz L002_I1.fastq.gz ...
                  --modality GEX
"""

import argparse
import gzip
import os
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# FASTQ helpers
# ---------------------------------------------------------------------------

def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def iter_records(fh):
    """Yield (header, seq, plus, qual) tuples."""
    while True:
        header = fh.readline()
        if not header:
            return
        yield header, fh.readline(), fh.readline(), fh.readline()


# ---------------------------------------------------------------------------
# Index detection
# ---------------------------------------------------------------------------

def scan_i1_multi(i1_paths, n_total=50000):
    """
    Pool I1 sequences from all lanes.
    Reads up to n_total/num_lanes reads from each lane so that every lane
    contributes equally to detection.
    """
    n_per_lane = max(1000, n_total // len(i1_paths))
    counter = Counter()
    for path in sorted(i1_paths):           # sort = consistent lane order
        with open_fastq(path) as fh:
            sampled = 0
            for _, seq, _, _ in iter_records(fh):
                counter[seq.strip()] += 1
                sampled += 1
                if sampled >= n_per_lane:
                    break
    return counter


def hamming(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))


def detect_samples(counter, run_id, min_freq=0.01):
    """
    Greedy centroid selection: any barcode above min_freq that is not within
    Hamming distance 1 of an already-chosen centroid becomes a new sample.
    Returns dict  centroid_seq -> sample_name  (descending frequency order).
    """
    total     = sum(counter.values())
    threshold = total * min_freq

    candidates = sorted(
        [(seq, cnt) for seq, cnt in counter.items() if cnt >= threshold],
        key=lambda x: -x[1],
    )

    centroids = {}
    for seq, _ in candidates:
        if any(hamming(seq, c) <= 1 for c in centroids):
            continue
        idx = len(centroids) + 1
        centroids[seq] = f"{run_id}_S{idx:02d}_{seq}"

    return centroids


def all_one_mismatch(seq):
    bases = "ACGTN"
    for i in range(len(seq)):
        for alt in bases:
            if alt != seq[i]:
                yield seq[:i] + alt + seq[i + 1:]


def build_lookup(centroids, mismatches=1):
    """
    barcode_seq -> sample_name dict including mismatch neighbours.
    Conflicting sequences (claimed by two samples) are dropped.
    """
    lookup    = {}
    conflicts = set()

    for centroid, sample_name in centroids.items():
        seqs = [centroid] + (list(all_one_mismatch(centroid)) if mismatches >= 1 else [])
        for seq in seqs:
            if seq in conflicts:
                continue
            if seq in lookup and lookup[seq] != sample_name:
                conflicts.add(seq)
                del lookup[seq]
            else:
                lookup[seq] = sample_name

    return lookup


# ---------------------------------------------------------------------------
# Demultiplexing
# ---------------------------------------------------------------------------

def demultiplex_lanes(lane_triples, lookup, sample_names, outdir):
    """
    Process every (r1, r2, i1) lane triple sequentially, writing all reads
    into the same per-sample output files.  Output files are opened once and
    kept open across all lanes, so memory usage stays flat.

    lane_triples : list of (r1_path, r2_path, i1_path) sorted by lane
    Returns      : Counter  sample_name -> read_count
    """
    os.makedirs(outdir, exist_ok=True)

    out_r1    = {s: gzip.open(os.path.join(outdir, f"{s}_S1_L001_R1_001.fastq.gz"), "wt")
                  for s in sample_names}
    out_r2    = {s: gzip.open(os.path.join(outdir, f"{s}_S1_L001_R2_001.fastq.gz"), "wt")
                  for s in sample_names}
    undet_r1  = gzip.open(os.path.join(outdir, "Undetermined_S0_L001_R1_001.fastq.gz"), "wt")
    undet_r2  = gzip.open(os.path.join(outdir, "Undetermined_S0_L001_R2_001.fastq.gz"), "wt")

    stats = Counter()

    for lane_idx, (r1_path, r2_path, i1_path) in enumerate(lane_triples, 1):
        print(f"[auto_demux]   Lane {lane_idx}/{len(lane_triples)}: {os.path.basename(r1_path)}",
              flush=True)
        with open_fastq(r1_path) as fh1, \
             open_fastq(r2_path) as fh2, \
             open_fastq(i1_path) as fi1:
            for r1_rec, r2_rec, i1_rec in zip(iter_records(fh1),
                                               iter_records(fh2),
                                               iter_records(fi1)):
                i1_seq = i1_rec[1].strip()
                sample = lookup.get(i1_seq)
                if sample:
                    out_r1[sample].writelines(r1_rec)
                    out_r2[sample].writelines(r2_rec)
                    stats[sample] += 1
                else:
                    undet_r1.writelines(r1_rec)
                    undet_r2.writelines(r2_rec)
                    stats["Undetermined"] += 1

    for fh in list(out_r1.values()) + list(out_r2.values()):
        fh.close()
    undet_r1.close()
    undet_r2.close()

    return stats


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_manifest(centroids, modality, outdir):
    path = os.path.join(outdir, "detected_samples.csv")
    with open(path, "w") as fh:
        fh.write("sample,centroid,modality\n")
        for centroid, name in centroids.items():
            fh.write(f"{name},{centroid},{modality}\n")


def write_stats(stats, outdir):
    total = sum(stats.values())
    path  = os.path.join(outdir, "demux_stats.csv")
    with open(path, "w") as fh:
        fh.write("sample,reads,fraction\n")
        for sample, count in sorted(stats.items(), key=lambda x: -x[1]):
            frac = count / total if total else 0.0
            fh.write(f"{sample},{count},{frac:.4f}\n")
            print(f"  {sample}: {count:,} reads  ({frac:.1%})", flush=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id",     required=True)
    ap.add_argument("--r1",         required=True, nargs="+",
                    help="R1 FASTQ(s) — one per lane, same order as --r2 / --i1")
    ap.add_argument("--r2",         required=True, nargs="+",
                    help="R2 FASTQ(s) — one per lane")
    ap.add_argument("--i1",         required=True, nargs="+",
                    help="I1 index FASTQ(s) — one per lane")
    ap.add_argument("--i2",         default=None, nargs="*",
                    help="I2 index FASTQ(s) — optional, one per lane")
    ap.add_argument("--modality",   required=True, help="GEX or TCR")
    ap.add_argument("--n-detect",   type=int,   default=50000)
    ap.add_argument("--min-freq",   type=float, default=0.01)
    ap.add_argument("--mismatches", type=int,   default=1)
    ap.add_argument("--outdir",     default="demuxed")
    args = ap.parse_args()

    # Sort all lane files by filename so lane order is always consistent
    r1_files = sorted(args.r1)
    r2_files = sorted(args.r2)
    i1_files = sorted(args.i1)

    n_lanes = len(r1_files)
    if not (len(r2_files) == len(i1_files) == n_lanes):
        sys.exit(f"[auto_demux] ERROR: --r1, --r2, --i1 must have the same number of files "
                 f"(got {n_lanes}, {len(r2_files)}, {len(i1_files)}).")

    lane_triples = list(zip(r1_files, r2_files, i1_files))

    # ------------------------------------------------------------------
    # Phase 1: detect barcodes (pooled across all lanes)
    # ------------------------------------------------------------------
    print(f"[auto_demux] Detecting barcodes from {n_lanes} lane(s), "
          f"sampling up to {args.n_detect:,} I1 reads total...", flush=True)
    counter   = scan_i1_multi(i1_files, n_total=args.n_detect)
    centroids = detect_samples(counter, args.run_id, min_freq=args.min_freq)

    if not centroids:
        sys.exit(
            f"[auto_demux] ERROR: No barcodes detected above "
            f"{args.min_freq * 100:.1f}% in {sum(counter.values()):,} sampled reads.\n"
            f"  → Try lowering --min-freq or increasing --n-detect.\n"
            f"  → Verify that --i1 files contain index reads (not transcript reads)."
        )

    total_sampled = sum(counter.values())
    print(f"[auto_demux] Detected {len(centroids)} sample(s) "
          f"from {total_sampled:,} sampled reads:", flush=True)
    for centroid, name in centroids.items():
        freq = counter[centroid] / total_sampled * 100
        print(f"  {name}  centroid={centroid}  freq={freq:.1f}%", flush=True)

    # ------------------------------------------------------------------
    # Phase 2: build lookup
    # ------------------------------------------------------------------
    lookup       = build_lookup(centroids, mismatches=args.mismatches)
    sample_names = list(centroids.values())

    # ------------------------------------------------------------------
    # Phase 3: manifest
    # ------------------------------------------------------------------
    os.makedirs(args.outdir, exist_ok=True)
    write_manifest(centroids, args.modality, args.outdir)

    # ------------------------------------------------------------------
    # Phase 4: demultiplex all lanes
    # ------------------------------------------------------------------
    print(f"[auto_demux] Demultiplexing {n_lanes} lane(s) into {args.outdir}/...",
          flush=True)
    stats = demultiplex_lanes(lane_triples, lookup, sample_names, args.outdir)

    # ------------------------------------------------------------------
    # Phase 5: stats
    # ------------------------------------------------------------------
    print("[auto_demux] Read assignment statistics:", flush=True)
    write_stats(stats, args.outdir)
    print("[auto_demux] Done.", flush=True)


if __name__ == "__main__":
    main()
