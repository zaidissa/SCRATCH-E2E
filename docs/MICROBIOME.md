# Microbiome stage — not implemented

`--run_microbiome` defaults to `false` and the pipeline errors with a clear
message if you force it on. This is not a merge decision; the stage was already
non-functional in the source repository.

## What is missing

`SCRATCH-QC-multi-mode/subworkflow/local/microbiome.nf` includes eight modules:

```
modules/local/microbiome/extract_unmapped/main.nf
modules/local/microbiome/bam_to_fastq/main.nf
modules/local/microbiome/filter_human_reads/main.nf
modules/local/microbiome/classify_microbiome/main.nf
modules/local/microbiome/refine_classification/main.nf
modules/local/microbiome/species_profiling/main.nf
modules/local/microbiome/link_microbes_to_cells/main.nf
modules/local/microbiome/generate_cell_profile/main.nf
```

None of them exist — `modules/local/microbiome/` is absent from the repository.
There was nothing to vendor.

Separately, the entry point `scratch_microbiome_entry.nf` includes
`./subworkflows/local/microbiome.nf` (plural), while the actual directory is
`subworkflow/` (singular), so the entry point could not have resolved either.

## The wiring that is in place

The DAG side is ready for these modules whenever they are written:

* `ALIGN` emits a `bam` channel — `(sample_id, possorted_genome_bam.bam)`
  extracted from the CellRanger `outs` directory. In the standalone pipeline the
  BAM was unreachable because the GEX and TCR alignment results shared a single
  channel variable.
* `main.nf` has a `microbiome` stage slot in the backbone/branch resolver, so
  `--from microbiome` and `--to microbiome` already parse.
* `params.input_bam_path` exists for entering the stage directly from BAMs.

## To implement it

1. Write the eight modules under `modules/local/microbiome/`.
2. Create `subworkflows/local/microbiome.nf` following the shape of the original
   (`EXTRACT_UNMAPPED → BAM_TO_FASTQ → FILTER_HUMAN_READS → CLASSIFY_MICROBIOME
   → REFINE_CLASSIFICATION`, with `SPECIES_PROFILING` and
   `LINK_MICROBES_TO_CELLS → GENERATE_CELL_PROFILE` in parallel).
3. Add `include { MICROBIOME } from './subworkflows/local/microbiome.nf'` and
   replace the `error` block in `main.nf` with the call, passing `ch_bam`.
4. Add container entries to `conf/containers.config` (KrakenUniq, Bracken,
   MetaPhlAn) and publish rules to `conf/modules.config`.
5. Add a `microbiome` branch to `INTEGRATE_OBJECTS` so per-cell microbial
   profiles join the master object like the other stages.
