# Numbat — allele-aware CNV calling

Opt-in fourth CNV caller. Off by default.

```bash
nextflow run . -profile docker \
    --run_numbat true \
    --numbat_gmap     /refs/numbat/genetic_map_hg38_withX.txt.gz \
    --numbat_snpvcf   /refs/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
    --numbat_paneldir /refs/numbat/1000G_hg38 \
    --input_samplesheet sheet.csv --input_exp_table meta.csv
```

## What it adds

inferCNV, SCEVAN and CopyKAT all infer copy number from **expression** smoothed
along the genome. That makes them structurally blind to **copy-neutral LOH**: an
event that changes allele ratios without changing dosage produces no expression
signal at all, so no amount of tuning will reveal it.

Numbat reads phased heterozygous SNPs out of the aligned BAM, so it sees the
allelic channel directly. It is not a better-tuned version of the other three —
it measures something they cannot.

Enabling it also makes `--stratify_min_callers 2` meaningful, since the
consensus then has three or four independent opinions instead of two.

## Requirements

| Requirement | Notes |
|---|---|
| Aligned BAM | `possorted_genome_bam.bam`, already emitted by `ALIGN` |
| Genetic map | Eagle2's `genetic_map_hg38_withX.txt.gz` |
| SNP VCF | e.g. `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf` |
| Phasing panel | 1000 Genomes reference panel, per-chromosome |
| Container | R `numbat` package + Eagle2 + samtools |

The three reference files are large external downloads and are **not** vendored
here. The Numbat project documents where to obtain them; fetch them once to a
shared location and point the three `--numbat_*` parameters at it. The pipeline
checks all three before starting and fails immediately with a named parameter if
one is missing — deliberately, because discovering a missing panel *after* a
pileup has spent hours walking a BAM is the worst possible time.

`--container_numbat` defaults to `pkgs/numbat:latest`. **Verify that tag, or
point it at your own build, before a production run** — it is the one image in
this pipeline I could not confirm.

## How it runs

```
ANNOTATION ──> NUMBAT_PREP ──> per-sample counts + barcode→cell_id map
                                        │
ALIGN (BAM) ────────────────────────────┴──> NUMBAT_PILEUP ──> phased allele counts
                                                                      │
                                             NUMBAT_RUN <─────────────┘
                                                  │
                                                  └──> *_numbat_calls.tsv ──> STRATIFY
```

Three processes rather than one, so the expensive pileup caches independently of
the model fit — retuning `--numbat_t` or `--numbat_max_iter` re-runs only
`NUMBAT_RUN`.

**The barcode problem, and why `NUMBAT_PREP` exists.** The BAM speaks raw
cellranger barcodes (`AACGT…-1`); the merged Seurat object speaks whatever
Seurat's merge produced (`SAMPLE_AACGT…-1`, `AACGT…-1_1`, …). Rather than
guessing that mapping with a regex at join time, `NUMBAT_PREP` writes an
explicit `barcode → cell_id` table per sample and `numbat_run.R` joins through
it, so the calls come back keyed to the pipeline's own cell ids. It also builds
the count matrix from **post-QC** cells, so Numbat does not model droplets that
were already discarded.

Samples are inner-joined on sample id: a sample needs both a BAM and surviving
post-QC cells to be processed. Samples that fail either condition are skipped
rather than failing the run.

## Parameters

| Parameter | Default | Notes |
|---|---|---|
| `--run_numbat` | `false` | |
| `--numbat_gmap` | — | required when enabled |
| `--numbat_snpvcf` | — | required when enabled |
| `--numbat_paneldir` | — | required when enabled |
| `--numbat_genome` | `hg38` | must match `--genome` |
| `--numbat_sample_col` | `sample_id` | column that defines a sample |
| `--numbat_min_cells` | `50` | below this a sample is skipped |
| `--numbat_t` | `1e-5` | transition probability |
| `--numbat_max_iter` | `2` | phylogeny refinement iterations |
| `--container_numbat` | `pkgs/numbat:latest` | verify before production use |

## How the calls are used

`numbat_run.R` prefers Numbat's own `compartment_opt` assignment and falls back
to thresholding `p_cnv` at 0.5. The result enters STRATIFY's evidence pool as a
fourth caller — no consensus logic changed, because the tally already counts
whatever callers ran.

The lineage gate still applies: a cell in a reference immune or stromal lineage
stays non-malignant even if Numbat calls it aneuploid. Numbat's disagreements
with the gate are counted in `data/cnv_false_positive_by_lineage.csv`, which is
a useful check on whether the allelic signal is behaving on your data.

## Caveats worth knowing before you commit

**Low-input samples are the real risk.** Numbat needs enough heterozygous-SNP
coverage per cell for the allelic channel to be informative. Low-cellularity,
low-tumour-content material — CSF aspirates in particular — is exactly where
allele-aware methods degrade. Expect it to behave well on resections and to be
marginal on CSF. Run it on one solid sample first and compare against
inferCNV/SCEVAN in the concordance table before extending to the whole cohort.

**The expression reference is generic.** `numbat_run.R` uses the package's
bundled `ref_hca` profile. A cohort-internal reference built from your own
immune cells would be preferable and is a worthwhile improvement, but it is not
available per-sample in the current wiring.

**On the benchmarking claim.** The mechanistic argument for allele-aware CNV
calling is solid and specific (copy-neutral LOH). My understanding is that
published benchmarks have generally favoured allele-aware methods, but I have
not verified the current literature — confirm it yourself before it drives a
methods decision.

**This path has been validated structurally, not scientifically.** The DAG,
channel joins, sample matching, config and failure modes are exercised by the
`-stub` run (50 tasks). The two R scripts have never executed against a real BAM
or a real Numbat installation. Treat the first real run as a shakedown: check
that `pileup_and_phase.R` accepts the argument names in your Numbat version, and
that `clone_post_*.tsv` carries `compartment_opt` or `p_cnv`.
