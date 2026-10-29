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
ANNOTATION ──> NUMBAT_PREP ──> per-INDIVIDUAL pooled counts + cell map
                                  │        + per-library barcodes
                                  │        + groups.tsv (library → individual)
                                  │
ALIGN (BAMs) ─────────────────────┴──> NUMBAT_PILEUP ──> phased allele counts
                                       one task per individual,      (one table
                                       all of their BAMs at once      per library)
                                                             │
                                       NUMBAT_RUN <──────────┘
                                       one task per individual
                                            │
                                            ├──> *_numbat_calls.tsv ──> STRATIFY
                                            ├──> *_numbat_by_library.tsv
                                            └──> numbat/<individual>/ ──> NUMBAT_REPORT
                                                                              │
                                                                              └──> figures
```

Three processes rather than one, so the expensive pileup caches independently of
the model fit — retuning `--numbat_t` or `--numbat_max_iter` re-runs only
`NUMBAT_RUN`.

## Numbat runs per individual, not per library

This is not a performance choice. **Numbat cannot report "there is no tumour
here."** `numbat:::mark_tumor_lineage()` assigns the tumour compartment with

```r
tumor_root  = mut_nodes[which.max(mut_burdens)]
compartment = ifelse(seq > 0, "tumor", "normal")
```

— no significance test, no null. Whatever segments survive `--numbat_min_llr`
go on a phylogeny and the heaviest clade becomes `tumor`, and `which.max()`
always returns something. Run per library, a PBMC or matched-normal library
therefore does **not** come back empty: it comes back with its own noise
partitioned into tumour and normal. That is exactly what was observed on
GBM-DFCI3 — neoplastic calls in libraries that should be ~0% tumour.

Pooled across the individual's libraries, the real aneuploid clones carry
overwhelming mutation burden, `tumor_root` lands on the genuine tumour clade,
and the normal cells stay in the normal compartment where they belong.

It is also what Numbat documents. Its own `pileup_and_phase.R`:

```
--label    Individual label. One per run.
--samples  Sample name(s); comma delimited if multiple.
           All samples must belong to the same individual.
```

Grouping therefore buys three things:

1. **A normal anchor**, as above.
2. **Joint phasing.** `numbat:::genotype(label, samples, …)` genotypes every
   library of the individual together, so a heterozygous site gets coverage from
   all of them — which is the allelic power Numbat has and the other callers
   do not.
3. **An internal expression reference.** `numbat::ref_hca` is 12 generic
   profiles from a different study, and the gap between it and the real cohort
   is read as expression CNV signal. With the libraries pooled, the
   individual's own immune and stromal cells are usually numerous enough to
   build a reference from instead. On a two-library test: 83 reference-lineage
   cells pooled (internal reference) versus 44 and 21 split per library (both
   fall back to `ref_hca`).

The grouping column is `--numbat_group_col`, default `patient_id`, falling back
through `patient`, `subject_id`, `donor_id` and finally the library column.
**That last fallback reproduces the per-library behaviour**, so it is reported as
a warning in both the task log and the run console rather than taken silently —
as is a cohort where every group happens to hold one library.

**Barcodes collide across libraries.** Two libraries of one patient both contain
`AACGT…-1`, so pooling them naively merges distinct cells. Every cell is
namespaced `<library>_<raw barcode>` — the same shape Numbat's own example uses
(`df_allele_example$cell` is `ATC2_TCCTAATCACAACCGC-1`). The pileup still writes
one allele table per library with raw barcodes; `numbat_run.R` applies the
prefix as it concatenates them.

**Why `NUMBAT_PREP` exists at all.** The BAM speaks raw cellranger barcodes; the
merged Seurat object speaks whatever Seurat's merge produced (`SAMPLE_AACGT…-1`,
`AACGT…-1_1`, …). Rather than guessing that mapping with a regex at join time,
`NUMBAT_PREP` writes an explicit `cell → barcode → cell_id` table and
`numbat_run.R` joins through it, so the calls come back keyed to the pipeline's
own cell ids. It also builds the count matrix from **post-QC** cells, so Numbat
does not model droplets that were already discarded.

`NUMBAT_PREP` runs in `image_cnv` (Seurat, no numbat) and `NUMBAT_RUN` in the
numbat image (numbat, no Seurat). That is why PREP *selects* the reference cells
and RUN *aggregates* them — `aggregate_counts()` needs the numbat package.

Libraries are inner-joined on library id: a library needs both a BAM and
surviving post-QC cells to be processed. Libraries that fail either condition
are skipped rather than failing the run.

### Checking the result

`NUMBAT_RUN` publishes `<individual>_numbat_by_library.tsv`:

| sample | n_cells | n_reference_lineage | malignant | pct_malignant_raw | malignant_gated | pct_malignant_gated | gate_overrode | reference_basis | n_reference_columns |
|---|---|---|---|---|---|---|---|---|---|
| …-TUMOUR-LIB | 4182 | 578 | 3604 | 86.2 | 3598 | 86.0 | 6 | internal | 11 |
| …-PBMC-LIB | 2910 | 2847 | 331 | 11.4 | 9 | 0.3 | 322 | internal | 11 |

**Read both percentages.** They answer different questions and only one of them
is the result.

`pct_malignant_raw` is Numbat's `compartment_opt` with nothing applied. For a
library with no malignant cells it is **structurally guaranteed to be non-zero**:
`mark_tumor_lineage()` always labels its heaviest clade `tumor`, so a PBMC
library returns its own noise partitioned into tumour and normal. A high number
in this column is not evidence of tumour.

`pct_malignant_gated` is what survives the lineage gate, and what the stratified
compartments actually use — a cell on one of the `--stratify_reference_labels`
lineages is forced non-malignant regardless of the CNV call. `gate_overrode` is
the difference, and reading it as an empirical false-positive rate for the CNV
calls is reasonable.

A PBMC library showing `11.4` raw and `0.3` gated is working correctly. A PBMC
library showing a high number in *both* is a real problem.

## What the expression reference is built from

`lambdas_ref` is a gene × pseudobulk matrix, and `run_numbat()` has **no
argument for reference cells** — no `ref_group_names`, no `normal_cells`. Unlike
inferCNV, you do not tell Numbat which cells are normal; you hand it a profile
matrix and it picks the best-fitting **column** per pseudobulk. Two consequences
follow directly:

- The cells that built the reference are **not held out**. Nothing connects
  `lambdas_ref` back to `count_mat`, so those same cells are still scored and can
  still be called tumour. The pipeline keeps them in `count_mat` deliberately —
  removing them would also remove the diploid anchor the phylogeny needs — and
  reports the gated percentage instead.
- A column's *name* is for us, not for Numbat.

### Columns are clusters, not cell type labels

The package author's own reference construction uses raw unsupervised clusters
with no annotation at all (`vignettes/mouse.Rmd` §4):

```r
clusters  = p2$clusters$PCA$multilevel
ref_annot = data.frame(cell = names(clusters), group = unname(clusters))
ref       = aggregate_counts(count_mat_ref,
                ref_annot %>% group_by(group) %>% filter(n() > 100))
```

What a column needs is transcriptional homogeneity and enough cells. scType
assigns **one label per cluster**, so grouping by label collapses several
clusters into one averaged column and discards resolution the clustering already
produced. `numbat_prep.R` therefore groups by `<label>__c<cluster>`, keeping the
same eligibility rule and giving Numbat finer columns to match against. Measured
on a fixture: 3 label columns → 4 cluster columns over the same 920 cells, with
one `TAM_Microglial` label contributing two distinct columns.

This also makes the per-column floor mean something. At label granularity every
column has hundreds of cells and a floor never fires; at cluster granularity it
is the actual quality control, which is why `--numbat_ref_min_group_cells` is
**100** — the author's number above — rather than 20.

### The five rules this follows

From Numbat's maintainer in [kharchenkolab/numbat#224][i224]:

| Rule | How it is satisfied |
|---|---|
| gene × cell-type matrix via `aggregate_counts()` **from raw counts**, *"don't use Seurat normalization"* | `numbat_prep.R` reads `LayerData(obj, layer = "counts")`; `numbat_run.R` calls `aggregate_counts(..., normalized = TRUE)` — raw in, normalized profile out |
| include the relevant normal cell types *present in your sample* | partially — see the caveat below |
| **no tumour/malignant cells** in the diploid reference | eligibility is `--stratify_reference_labels`, which excludes `Malignant` and the four labels Neftel malignant states also carry |
| internal normals beat external datasets | the reference is built from the individual's own cells; `numbat::ref_hca` is only a fallback, and the log says when it fired |
| raise `min_LLR` / lower `max_entropy` against false positives in normal cells | both reachable as `--numbat_min_llr` and `--numbat_max_entropy`, both left at Numbat's own defaults |

[i224]: https://github.com/kharchenkolab/numbat/issues/224

**The caveat on rule 2, stated plainly.** GBM tissue contains normal astrocytes,
OPCs and neurons, and the GBM reference has **no glial or neural column at all**.
That is not an oversight: the Neftel AC-like, OPC-like and NPC-like malignant
states carry exactly the `Astrocyte`, `OPC`, `Radial_Glia` and `Neuron` labels,
so admitting those labels would admit tumour and break rule 3. No label-based
rule can separate the two.

A *cluster* can be judged individually, with evidence from outside the
annotation — a flat inferCNV profile, say. `--numbat_extra_reference_clusters`
admits named clusters regardless of label, for exactly that. It is empty by
default and should stay that way unless you have that evidence: deciding it from
the annotation alone is circular, because the annotation is what the CNV call is
meant to adjudicate.

### Seeing what the baseline was

`NUMBAT_PREP` publishes `<individual>_refgroups.tsv`, one row per column,
including the ones the floor dropped:

| ref_group | celltype | cluster | n_cells | admitted | basis |
|---|---|---|---|---|---|
| Monocyte__c6 | Monocyte | 6 | 250 | yes | reference_lineage |
| TAM_Microglial__c1 | TAM_Microglial | 1 | 150 | yes | reference_lineage |
| T_Cells__c4 | T_Cells | 4 | 40 | no_below_floor | reference_lineage |

Without it the only trace of a dropped column is a line in `.command.log`.

## Parameters

| Parameter | Default | Notes |
|---|---|---|
| `--run_numbat` | `false` | |
| `--numbat_gmap` | — | required when enabled |
| `--numbat_snpvcf` | — | required when enabled |
| `--numbat_paneldir` | — | required when enabled |
| `--numbat_genome` | `hg38` | must match `--genome` |
| `--numbat_sample_col` | `sample_id` | column that defines a **library** (must match the BAM directory names) |
| `--numbat_group_col` | `patient_id` | column that defines an **individual** — the unit Numbat runs on |
| `--numbat_min_ref_cells` | `100` | reference-lineage cells needed before an internal expression reference replaces `ref_hca`. Raised from 50 together with the per-column floor: a column can never exceed the total |
| `--numbat_cluster_col` | `seurat_clusters` | cluster column the reference columns are built on. Unset/unmatched degrades to one column per cell type, with a warning |
| `--numbat_ref_min_group_cells` | `100` | per-**column** cell floor; the author's number in `mouse.Rmd` §4 |
| `--numbat_extra_reference_clusters` | `''` | clusters admitted to the reference regardless of label. Empty by default — see the caveat above |
| `--numbat_min_cells` | `50` | below this an individual is skipped |
| `--numbat_t` | `1e-5` | transition probability |
| `--numbat_max_iter` | `2` | phylogeny refinement iterations |
| `--numbat_min_llr` | `5` | log-likelihood-ratio floor for keeping a CNV. Numbat's own default; raising it is the maintainer's first lever against false positives in normal cells |
| `--numbat_max_entropy` | `0.5` | single-cell entropy ceiling. Numbat's own default; lowering it is the second lever. Neither has been calibrated on this cohort, so both are left where Numbat puts them |
| `--skip_numbat_figures` | `false` | turn off `NUMBAT_REPORT` |
| `--numbat_figures_ncut` | `4` | clone counts to re-cut the phylogeny at; `0` skips that panel |
| `--numbat_figures_min_clone_cells` | `50` | smaller clones are left out of the pseudobulk figure |
| `--numbat_figures_min_llr` | `10` | LLR floor for drawing a CNV on the pseudobulk figure |
| `--container_numbat` | `pkgs/numbat:latest` | verify before production use |

## Figures

`run_numbat()` is called with `plot = FALSE`, so Numbat itself draws nothing.
`NUMBAT_REPORT` draws the figures afterwards instead, per individual, from the
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
