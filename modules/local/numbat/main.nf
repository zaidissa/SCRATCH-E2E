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
        tuple val(sample_id), path(bam), path(barcodes)
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

        # Cellranger emits the index alongside the BAM; regenerate if absent.
        if [ ! -f "${bam}.bai" ]; then
            samtools index -@ ${task.cpus} ${bam}
        fi

        pileup_and_phase.R \\
            --label     ${sample_id} \\
            --samples   ${sample_id} \\
            --bams      ${bam} \\
            --barcodes  ${barcodes} \\
            --outdir    pileup \\
            --gmap      ${gmap} \\
            --snpvcf    ${snpvcf} \\
            --paneldir  ${paneldir} \\
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

        numbat_run.R \\
            ${counts} ${alleles} ${cellmap} \\
            ${sample_id} numbat ${task.cpus} \\
            ${params.numbat_min_cells} ${params.numbat_t} ${params.numbat_max_iter}

        # Numbat's own diagnostic plots, when it produced any.
        find numbat/${sample_id} -name "*.png" -exec cp {} figures/numbat/ \\; 2>/dev/null || true
        """

    stub:
        """
        mkdir -p numbat/${sample_id} figures/numbat
        printf 'barcode\\tcell_id\\tsample\\tnumbat_call\\tp_cnv\\n' > numbat/${sample_id}_numbat_calls.tsv
        printf 'AACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\tmalignant\\t0.98\\n' >> numbat/${sample_id}_numbat_calls.tsv
        touch numbat/${sample_id}/stub.txt
        """
}
