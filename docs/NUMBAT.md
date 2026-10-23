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

inferCNV and CopyKAT infer copy number from **expression** smoothed
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

### Where the reference files come from

Not vendored here — they are large external downloads. Fetch once to a shared
location (on HPC put them somewhere group-readable) and point the three
`--numbat_*` parameters at them. `bin/fetch_numbat_refs.sh` does all of this.

| File | Source | Size |
|---|---|---|
| **SNP VCF** `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz` | `https://sourceforge.net/projects/cellsnp/files/SNPlist/` — the cellsnp-lite SNP lists Numbat's docs point to. Swap `hg38`→`hg19` for the other build. | ~200 MB |
| **Phasing panel** `1000G_hg38/` | `http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip` — the Kharchenko lab's prepared 1000G panel (`1000G_hg19.zip` for hg19). Unzip; point `--numbat_paneldir` at the resulting directory. | ~9 GB unzipped |
| **Genetic map** `genetic_map_hg38_withX.txt.gz` | Ships **inside Eagle2**: `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz` → `Eagle_v2.4.1/tables/`. Also already present in the Numbat container at `/Eagle_v2.4.1/tables/`. | ~15 MB |
| **BAM** `possorted_genome_bam.bam` | Produced by this pipeline's own `ALIGN` stage, or your existing CellRanger `outs/`. Must be the position-sorted BAM **with its `.bai`**. | per sample |

All four URLs were checked and resolved at the time of writing; if one has moved,
the canonical index is the Numbat repository: `https://github.com/kharchenkolab/numbat`.

The genome build must be consistent across all of them and match `--numbat_genome`
and `--genome`. Mixing an hg19 panel with hg38 BAMs will not error cleanly — it
produces garbage calls, which is far worse.

The pipeline checks all three paths before starting and fails immediately naming
the missing parameter — deliberately, because discovering a missing panel *after*
a pileup has spent hours walking a BAM is the worst possible time.

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
                                                  ├──> *_numbat_calls.tsv ──> STRATIFY
                                                  │
                                                  └──> numbat/<sample>/ ──> NUMBAT_REPORT
                                                                                  │
                                                                                  └──> figures
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
| `--numbat_min_llr` | `5` | log-likelihood-ratio floor for keeping a CNV |
| `--skip_numbat_figures` | `false` | turn off `NUMBAT_REPORT` |
| `--numbat_figures_ncut` | `4` | clone counts to re-cut the phylogeny at; `0` skips that panel |
| `--numbat_figures_min_clone_cells` | `50` | smaller clones are left out of the pseudobulk figure |
| `--numbat_figures_min_llr` | `10` | LLR floor for drawing a CNV on the pseudobulk figure |
| `--container_numbat` | `pkgs/numbat:latest` | verify before production use |

## Figures

`run_numbat()` is called with `plot = FALSE`, so Numbat itself draws nothing.
`NUMBAT_REPORT` draws the figures afterwards instead, per sample, from the
output directory `NUMBAT_RUN` has already written — the posteriors, consensus
segments, phylogeny and mutation graph are all written outside `run_numbat`'s
own `if (plot)` branch, so nothing had to be recomputed to get them.

It is a separate process on purpose. Nextflow's task hash covers a process's
script text, so adding one plotting line to `NUMBAT_RUN` re-runs the model fit
for every sample — 957 CPU-hours over 15 GBM libraries. Turning `plot = TRUE`
on would do the same, *and* produce only five of the figures below: Numbat's
internal plotting draws no phylogeny heatmap, no single-cell tree, no mutation
history and no consensus plot.

Per sample, published to `cnv/figures/` as PNG (plus PDF where the figure is
not too dense to vectorise):

| File | What it shows |
|---|---|
| `<sample>_phylo_heatmap.png` | the signature figure: single-cell phylogeny beside the per-cell CNV calls, with a clone bar and the tumour/normal divide |
| `<sample>_phylo_heatmap_by_ncut.png` | the same heatmap re-cut at 1–4 clones, because the clone count is a choice rather than a measurement |
| `<sample>_consensus_segments.png` | the CNV events along the genome |
| `<sample>_segment_heatmap.png` | the same segments as a genome-wide heatmap, one row per group |
| `<sample>_bulk_clone_profiles.png` | pseudobulk expression log-ratio over phased allele frequency, per clone — the only figure that shows the allelic evidence directly |
| `<sample>_sc_phylogeny.png` | the phylogeny alone, each CNV placed on the branch where it was acquired |
| `<sample>_mutation_history.png` | the clone tree, nodes sized by clone |
| `<sample>_expression_rolling.png` | window-smoothed expression, clustered — the expression channel on its own |

`cnv/data/<sample>_numbat_figure_manifest.tsv` names every figure, whether it
was drawn, and why not when it was not. A figure browser cannot tell "not
drawn" from "nothing to draw"; the manifest can. Every figure is attempted
independently, and `NUMBAT_REPORT` is `errorStrategy = 'ignore'` — a picture
must never cost a call.

The CNV concordance report adds three more, cohort-wide, on the pipeline's own
UMAP: Numbat's clone assignment; its aneuploidy posterior split into the
**expression** and **allele** channels beside the joint result; and the most
widely carried CNV segments one panel at a time. The split is the diagnostic —
a cell that is positive on the allele channel but not on expression carries a
copy-neutral event that no expression-based caller can see at all.

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
inferCNV in the concordance table before extending to the whole cohort.

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
