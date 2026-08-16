/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Numbat — allele-aware CNV calling  (opt-in)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Enabled with --run_numbat true. Off by default because, unlike the other
    three callers, it needs the aligned BAM plus a phasing reference panel.

    What it adds over inferCNV / SCEVAN / CopyKAT: those infer copy number from
    expression smoothed along the genome, so they are structurally blind to
    copy-neutral LOH — an event that changes allele ratios without changing
    dosage produces no expression signal. Numbat reads phased heterozygous SNPs
    from the BAM and sees it.

    Three processes, kept separate so the expensive pileup caches independently
    of the model fit:

      NUMBAT_PREP     annotated object -> per-sample counts + an explicit
                      barcode -> cell_id map (the BAM speaks raw cellranger
                      barcodes; the merged object does not)
      NUMBAT_PILEUP   BAM + barcodes  -> phased allele counts   [per sample]
      NUMBAT_RUN      counts + alleles -> per-cell calls         [per sample]

    The calls join STRATIFY's evidence pool as a fourth caller. No consensus
    logic changes: the tally already counts whatever callers ran, so enabling
    Numbat is also what makes --stratify_min_callers 2 meaningful.
----------------------------------------------------------------------------------------
*/

// Fetches the three external references Numbat needs (~10 GB) when they were
// not supplied. `storeDir` is what makes this safe to leave on: the files live
// OUTSIDE work/, so they survive -resume, a cleaned work directory, and every
// later run -- the download happens exactly once per destination, and the
// process is skipped entirely whenever the outputs are already there.
//
// The build is part of the stored path on purpose. An hg19 panel against hg38
// BAMs does not error, it silently produces wrong calls, so hg19 and hg38 must
// never be able to occupy the same cache slot.
process NUMBAT_FETCH_REFS {

    tag "Numbat refs (${params.numbat_genome})"
    label 'process_low'
    storeDir "${params.numbat_refs_dir}/${params.numbat_genome}"

    output:
        path "genetic_map_${params.numbat_genome}_withX.txt.gz", emit: gmap
        path "genome1K.phase3.SNP_AF5e2.chr1toX.${params.numbat_genome}.vcf.gz", emit: snpvcf
        path "1000G_${params.numbat_genome}", type: 'dir', emit: paneldir

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        fetch_numbat_refs.sh . ${params.numbat_genome}
        # The fetch script is idempotent and verifies each file exists, but a
        # truncated download would still satisfy `-f`. Fail here rather than let
        # a partial panel produce quietly wrong phasing.
        for f in genetic_map_${params.numbat_genome}_withX.txt.gz \
                 genome1K.phase3.SNP_AF5e2.chr1toX.${params.numbat_genome}.vcf.gz; do
            [ -s "\$f" ] || { echo "numbat reference \$f is missing or empty" >&2; exit 1; }
        done
        [ -d "1000G_${params.numbat_genome}" ] || {
            echo "numbat phasing panel directory is missing" >&2; exit 1; }
        """

    stub:
        """
        touch genetic_map_${params.numbat_genome}_withX.txt.gz
        touch genome1K.phase3.SNP_AF5e2.chr1toX.${params.numbat_genome}.vcf.gz
        mkdir -p 1000G_${params.numbat_genome}
        """
}


process NUMBAT_PREP {

    tag "Numbat: per-sample inputs"
    label 'process_medium'

    input:
        path seurat_object

    output:
        path "numbat_input/*_counts.rds",   emit: counts
        path "numbat_input/*_barcodes.tsv", emit: barcodes
        path "numbat_input/*_cellmap.tsv",  emit: cellmap
        path "numbat_input/samples.txt",    emit: samples

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p numbat_input
        numbat_prep.R ${seurat_object} ${params.numbat_sample_col} numbat_input
        """

    stub:
        """
        mkdir -p numbat_input
        touch numbat_input/TESTA_counts.rds
        echo "AACGTAACGTAACGTA-1" > numbat_input/TESTA_barcodes.tsv
        printf 'cell_id\\tbarcode\\tsample\\nTESTA_AACGTAACGTAACGTA-1\\tAACGTAACGTAACGTA-1\\tTESTA\\n' > numbat_input/TESTA_cellmap.tsv
        echo TESTA > numbat_input/samples.txt
        """
}


process NUMBAT_PILEUP {

    tag "Numbat pileup: ${sample_id}"
    label 'process_high'

    // The heavy step: walks the BAM for allele counts at heterozygous SNPs and
    // phases them against a population reference panel.

    input:
        tuple val(sample_id), path(bam), path(bai), path(barcodes)
        path  gmap
        path  snpvcf
        path  paneldir

    output:
        tuple val(sample_id), path("pileup/${sample_id}_allele_counts.tsv.gz"), emit: alleles

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p pileup

        # Reference resolution. The numbat image already ships all three (verified
        # in pkharchenkolab/numbat-rbase: /data/1000G_<build>, /data/genome1K...vcf,
        # /Eagle_v2.4.1/tables/genetic_map_<build>_withX.txt.gz), so the default is
        # to use those. Downloading them is a fallback, not the happy path:
        # the upstream panel URL (pklab.med.harvard.edu) sits behind an Incapsula
        # WAF that answers automated requests with a 212-byte HTML challenge and
        # HTTP 200, so curl saves the challenge page AS the zip and unzip fails.
        #
        # A staged NO_FILE placeholder means "not supplied -- use the image copy".
        GMAP="${gmap}"; SNPVCF="${snpvcf}"; PANELDIR="${paneldir}"
        case "\$GMAP"     in *NO_FILE*) GMAP=/Eagle_v2.4.1/tables/genetic_map_${params.numbat_genome}_withX.txt.gz ;; esac
        case "\$SNPVCF"   in *NO_FILE*) SNPVCF=\$(ls /data/genome1K.phase3.SNP_AF5e2.chr1toX.${params.numbat_genome}.vcf* 2>/dev/null | head -1) ;; esac
        case "\$PANELDIR" in *NO_FILE*) PANELDIR=/data/1000G_${params.numbat_genome} ;; esac

        # storeDir caches on EXISTENCE, not validity, and a WAF challenge page or
        # a truncated download both satisfy `-f`. Check before doing any work.
        for ref in "\$GMAP" "\$SNPVCF"; do
            [ -s "\$ref" ] || { echo "numbat reference '\$ref' is missing or empty." >&2; exit 1; }
        done
        if [ -z "\$(ls -A "\$PANELDIR" 2>/dev/null)" ]; then
            echo "numbat phasing panel '\$PANELDIR' is empty or missing." >&2; exit 1
        fi
        echo "[numbat] gmap=\$GMAP"; echo "[numbat] snpvcf=\$SNPVCF"; echo "[numbat] paneldir=\$PANELDIR"

        # Cellranger emits the index alongside the BAM; regenerate if absent.
        # pileup_and_phase.R is NOT on PATH in pkharchenkolab/numbat-rbase -- it
        # ships inside the installed R package at
        # /usr/local/lib/R/site-library/numbat/bin/. Ask R where it is rather than
        # hardcoding the path, so a different image or version still resolves.
        PILEUP=\$(command -v pileup_and_phase.R || true)
        if [ -z "\$PILEUP" ]; then
            PILEUP=\$(Rscript -e 'cat(system.file("bin/pileup_and_phase.R", package = "numbat"))' 2>/dev/null || true)
        fi
        [ -n "\$PILEUP" ] && [ -s "\$PILEUP" ] || {
            echo "pileup_and_phase.R not found on PATH or in the numbat R package." >&2; exit 1; }
        echo "[numbat] pileup_and_phase.R=\$PILEUP"

        # Re-indexing an 80 GB BAM costs more than the pileup, so the .bai is
        # staged alongside it; only build one if it genuinely is not there.
        if [ ! -s "${bam}.bai" ]; then
            echo "[numbat] no .bai staged -- indexing ${bam} (slow)" >&2
            samtools index -@ ${task.cpus} ${bam}
        fi

        Rscript \$PILEUP \\
            --label     ${sample_id} \\
            --samples   ${sample_id} \\
            --bams      ${bam} \\
            --barcodes  ${barcodes} \\
            --outdir    pileup \\
            --gmap      \$GMAP \\
            --snpvcf    \$SNPVCF \\
            --paneldir  \$PANELDIR \\
            --ncores    ${task.cpus}

        # pileup_and_phase.R names its output <label>_allele_counts.tsv.gz;
        # normalise in case a version writes it elsewhere in the tree.
        if [ ! -f "pileup/${sample_id}_allele_counts.tsv.gz" ]; then
            found=\$(find pileup -name "*allele_counts.tsv.gz" | head -1)
            [ -n "\$found" ] && cp "\$found" "pileup/${sample_id}_allele_counts.tsv.gz"
        fi
        """

    stub:
        """
        mkdir -p pileup
        printf 'cell\\tsnp_id\\tCHROM\\tPOS\\tREF\\tALT\\tAD\\tDP\\tGT\\n' | gzip > pileup/${sample_id}_allele_counts.tsv.gz
        """
}


process NUMBAT_RUN {

    tag "Numbat: ${sample_id}"
    label 'process_high'

    input:
        tuple val(sample_id), path(alleles), path(counts), path(cellmap)

    output:
        path "numbat/${sample_id}_numbat_calls.tsv", emit: calls
        // The whole point of the nocalls note is to explain a 0-row table to
        // whoever reads the RESULTS, so it has to be a declared output --
        // undeclared files stay in the work dir and publishDir copies nothing.
        path "numbat/${sample_id}_numbat_nocalls.txt", emit: nocalls, optional: true
        path "numbat/${sample_id}/**",               emit: detail,  optional: true
        path "figures/numbat/**",                    emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p numbat figures/numbat

        export NUMBAT_GENOME=${params.numbat_genome}
        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}

        # Echo the threshold BEFORE running. Nextflow does not hash bin/ script
        # contents, so an edit there alone will not re-run a cached task -- an
        # earlier 3-hour run silently reused a stale task and reported a result
        # for a threshold it never used. This line makes the value visible in
        # .command.log within seconds, and its presence in the script also
        # changes the task hash so the fix actually takes effect.
        echo "[numbat] min_LLR=${params.numbat_min_llr} min_cells=${params.numbat_min_cells} t=${params.numbat_t} max_iter=${params.numbat_max_iter}"

        numbat_run.R \\
            ${counts} ${alleles} ${cellmap} \\
            ${sample_id} numbat ${task.cpus} \\
            ${params.numbat_min_cells} ${params.numbat_t} ${params.numbat_max_iter} \\
            ${params.numbat_min_llr}

        # Numbat's own diagnostic plots, when it produced any.
        find numbat/${sample_id} -name "*.png" -exec cp {} figures/numbat/ \\; 2>/dev/null || true
        """

    stub:
        """
        mkdir -p numbat/${sample_id} figures/numbat
        printf 'barcode\\tcell_id\\tsample\\tnumbat_call\\tp_cnv\\n' > numbat/${sample_id}_numbat_calls.tsv
        printf 'AACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\tmalignant\\t0.98\\n' >> numbat/${sample_id}_numbat_calls.tsv
        printf 'stub: null-result note\\n' > numbat/${sample_id}_numbat_nocalls.txt
        touch numbat/${sample_id}/stub.txt
        """
}
