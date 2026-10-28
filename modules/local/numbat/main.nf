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

    tag "Numbat: per-individual inputs"
    label 'process_medium'

    // Splits the annotated object by INDIVIDUAL, not by library. See the header
    // of bin/numbat_prep.R for why: numbat:::mark_tumor_lineage() has no null
    // and always labels its heaviest clade "tumor", so a library with no
    // malignant cells cannot come back empty -- it comes back with its own
    // noise partitioned into tumour and normal. Pooling an individual's
    // libraries gives the real clones something to dominate.
    //
    // `helpers` is assets/viz/scratch_report.R, staged for its label matcher.
    // The reference lineages that define the null in STRATIFY must be the same
    // cells that define the expression reference here, and two copies of one
    // matching rule diverge silently -- the notebook that draws the concordance
    // figure already shares it for exactly that reason.

    input:
        path seurat_object
        path helpers

    output:
        path "numbat_input/*_counts.rds",        emit: counts
        path "numbat_input/*_barcodes.tsv",      emit: barcodes
        path "numbat_input/*_cellmap.tsv",       emit: cellmap
        path "numbat_input/*_refcells.tsv",      emit: refcells
        path "numbat_input/groups.tsv",          emit: groups
        path "numbat_input/grouping_report.tsv", emit: grouping_report
        path "numbat_input/samples.txt",         emit: samples

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p numbat_input
        numbat_prep.R ${seurat_object} \\
            ${params.numbat_sample_col} \\
            ${params.numbat_group_col} \\
            ${params.celltype_col} \\
            '${params.stratify_reference_labels}' \\
            ${params.numbat_min_ref_cells} \\
            numbat_input ${helpers}
        """

    stub:
        """
        mkdir -p numbat_input
        touch numbat_input/PT1_counts.rds
        printf 'cell\\tcelltype\\n' > numbat_input/PT1_refcells.tsv
        echo "AACGTAACGTAACGTA-1" > numbat_input/TESTA_barcodes.tsv
        echo "AACGTAACGTAACGTC-1" > numbat_input/TESTB_barcodes.tsv
        printf 'cell_id\\tbarcode\\tcell\\tlibrary\\tgroup\\n' > numbat_input/PT1_cellmap.tsv
        printf 'TESTA_AACGTAACGTAACGTA-1\\tAACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\tPT1\\n' >> numbat_input/PT1_cellmap.tsv
        printf 'TESTB_AACGTAACGTAACGTC-1\\tAACGTAACGTAACGTC-1\\tTESTB_AACGTAACGTAACGTC-1\\tTESTB\\tPT1\\n' >> numbat_input/PT1_cellmap.tsv
        printf 'library\\tgroup\\nTESTA\\tPT1\\nTESTB\\tPT1\\n' > numbat_input/groups.tsv
        printf 'group\\tn_libraries\\tlibraries\\tn_cells\\tn_reference_cells\\treference_basis\\nPT1\\t2\\tTESTA,TESTB\\t2\\t0\\tref_hca\\n' > numbat_input/grouping_report.tsv
        echo PT1 > numbat_input/samples.txt
        """
}


process NUMBAT_PILEUP {

    tag "Numbat pileup: ${group_id}"
    label 'process_high'

    // The heavy step: walks the BAM for allele counts at heterozygous SNPs and
    // phases them against a population reference panel.

    // ONE TASK PER INDIVIDUAL, not per library. pileup_and_phase.R takes
    // comma-delimited --samples/--bams/--barcodes and says so itself: "All
    // samples must belong to the same individual." It genotypes them jointly
    // (numbat:::genotype(label, samples, ...)), so a heterozygous site gets
    // coverage from every library the individual contributed -- which is the
    // allelic power that expression-only callers have none of.
    //
    // stageAs is REQUIRED here and not cosmetic: every cellranger BAM is named
    // possorted_genome_bam.bam, so staging several into one task directory
    // collides and Nextflow aborts with "input file name collision". 'bam??/*'
    // puts each in its own numbered directory keeping its basename, and the
    // .bai lands beside its .bam because both declarations number in list
    // order. Measured on 25.04.8 (Cirro's version) and 26.04.6: three BAMs
    // passed in the order C,A,B staged as in01=C, in02=A, in03=B, with each
    // index in the matching directory -- list order, not sorted. The subworkflow
    // sorts the rows by library before building the lists, so the order is
    // deterministic and the task hash is stable across runs.

    input:
        tuple val(group_id), val(libs),
              path(bams,     stageAs: 'bam??/*'),
              path(bais,     stageAs: 'bam??/*'),
              path(barcodes, stageAs: 'bc??/*')
        path  gmap
        path  snpvcf
        path  paneldir

    output:
        // One allele table PER LIBRARY, jointly phased. numbat_run.R
        // concatenates them, namespacing each cell as "<library>_<barcode>".
        tuple val(group_id), path("pileup/*_allele_counts.tsv.gz"), emit: alleles

    when:
        task.ext.when == null || task.ext.when

    script:
        def libList = (libs      instanceof List ? libs      : [libs]).collect { it.toString() }
        def bamList = (bams      instanceof List ? bams      : [bams]).collect { it.toString() }
        def bcList  = (barcodes  instanceof List ? barcodes  : [barcodes]).collect { it.toString() }
        if (libList.size() != bamList.size() || libList.size() != bcList.size())
            error "NUMBAT_PILEUP ${group_id}: ${libList.size()} librar(ies), " +
                  "${bamList.size()} BAM(s), ${bcList.size()} barcode file(s) -- these must match."
        """
        set -euo pipefail
        mkdir -p pileup

        echo "[numbat] individual=${group_id} libraries=${libList.join(',')}"

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
        #
        # LD_LIBRARY_PATH, scoped to this one command. The image carries TWO
        # htslibs: Debian's 1.20 under /lib/x86_64-linux-gnu, which /usr/bin/
        # samtools 1.20 is linked against, and a hand-built 1.15.1 under
        # /usr/local/lib, which cellsnp-lite is linked against. /usr/local/lib
        # wins in ld.so.conf, so samtools loads the 1.15.1 copy and dies with
        #     samtools: symbol lookup error: undefined symbol: bam_aux_first,
        #     version HTSLIB_1.17
        # Upstream pkharchenkolab/numbat-rbase does NOT have this problem even
        # though both libraries are present and the samtools binary is byte
        # identical: its ld.so.cache was built before /usr/local/lib held an
        # htslib. Adding procps in containers/numbat/Dockerfile ran apt's
        # ldconfig trigger, the cache was regenerated, and /usr/local/lib took
        # over -- so the procps layer broke samtools without touching it.
        # Measured by diffing the two images (Cirro 2026-10-21, exit 127).
        #
        # Deleting the old library is NOT the fix: cellsnp-lite, which does the
        # actual pileup, resolves to /usr/local/lib/libhts.so.3. Hence the
        # per-command override rather than an ENV or a global ldconfig.
        for b in ${bamList.join(' ')}; do
            if [ ! -s "\$b.bai" ]; then
                echo "[numbat] no .bai staged for \$b -- indexing (slow)" >&2
                LD_LIBRARY_PATH=/lib/x86_64-linux-gnu samtools index -@ ${task.cpus} "\$b"
            fi
        done

        Rscript \$PILEUP \\
            --label     ${group_id} \\
            --samples   ${libList.join(',')} \\
            --bams      ${bamList.join(',')} \\
            --barcodes  ${bcList.join(',')} \\
            --outdir    pileup \\
            --gmap      \$GMAP \\
            --snpvcf    \$SNPVCF \\
            --paneldir  \$PANELDIR \\
            --ncores    ${task.cpus}

        # pileup_and_phase.R writes {outdir}/{sample}_allele_counts.tsv.gz once
        # per library. Fail loudly if a library produced none: an allele table
        # that is quietly absent costs that library its allelic channel, which is
        # the only thing Numbat has over inferCNV, and nothing downstream can
        # tell the difference between "no table" and "no events".
        for lb in ${libList.join(' ')}; do
            if [ ! -s "pileup/\${lb}_allele_counts.tsv.gz" ]; then
                found=\$(find pileup -name "\${lb}*allele_counts.tsv.gz" | head -1)
                if [ -n "\$found" ]; then
                    cp "\$found" "pileup/\${lb}_allele_counts.tsv.gz"
                else
                    echo "pileup produced no allele table for library \$lb of ${group_id}." >&2
                    exit 1
                fi
            fi
        done
        """

    stub:
        // `script:` and `stub:` are separate closures -- a def in one is not
        // visible in the other ("No such variable: libList"), and because a stub
        // only runs under -stub the failure hides until someone stub-tests the DAG.
        def libList = (libs instanceof List ? libs : [libs]).collect { it.toString() }
        """
        mkdir -p pileup
        for lb in ${libList.join(' ')}; do
            printf 'cell\\tsnp_id\\tCHROM\\tPOS\\tREF\\tALT\\tAD\\tDP\\tGT\\n' | gzip > "pileup/\${lb}_allele_counts.tsv.gz"
        done
        """
}


process NUMBAT_RUN {

    tag "Numbat: ${group_id}"
    label 'process_high'

    // One task per INDIVIDUAL, over that individual's pooled libraries. The
    // allele tables arrive one per library (jointly phased) and keep their
    // library-named basenames, which is how numbat_run.R knows which prefix to
    // put on each cell before concatenating them.

    input:
        tuple val(group_id), path(alleles), path(counts), path(cellmap), path(refcells)

    output:
        path "numbat/${group_id}_numbat_calls.tsv", emit: calls
        // The whole point of the nocalls note is to explain a 0-row table to
        // whoever reads the RESULTS, so it has to be a declared output --
        // undeclared files stay in the work dir and publishDir copies nothing.
        path "numbat/${group_id}_numbat_nocalls.txt", emit: nocalls, optional: true
        // Malignant fraction PER LIBRARY within the individual. This is the
        // number the old per-library configuration got wrong and the first one a
        // reviewer checks: a PBMC or matched-normal library should be ~0%
        // neoplastic, and now that is a published table rather than something
        // only visible after plotting the final object.
        path "numbat/${group_id}_numbat_by_library.tsv", emit: by_library, optional: true
        path "numbat/${group_id}/**",               emit: detail,  optional: true
        // The EVENT-level table: one row per CNV segment, with cnv_state
        // (loh/amp/del/bamp), coordinates, LLR and gene/SNP support. The
        // per-cell calls above answer "how many cells", not "how many CNVs" --
        // a distinction worth keeping, since they differ by two orders of
        // magnitude (9,227 cells vs 66 segments on this data).
        path "numbat/${group_id}/segs_consensus_*.tsv", emit: segments, optional: true
        path "figures/numbat/**",                   emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def alleleList = (alleles instanceof List ? alleles : [alleles]).collect { it.toString() }
        """
        set -euo pipefail
        mkdir -p numbat figures/numbat

        export NUMBAT_GENOME=${params.numbat_genome}
        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}

        # Echo the parameters BEFORE running. bin/ IS hashed per referenced
        # script -- editing bin/numbat_run.R does re-run this task, measured on
        # 25.04.8 and 26.04.6 -- but a params change alone would not be visible
        # anywhere, and an earlier 3-hour run reported a result for a threshold
        # it never used. This line puts the values in .command.log within
        # seconds of the task starting.
        echo "[numbat] group=${group_id} libraries=${alleleList.size()} min_LLR=${params.numbat_min_llr} min_cells=${params.numbat_min_cells} t=${params.numbat_t} max_iter=${params.numbat_max_iter}"

        numbat_run.R \\
            ${counts} ${cellmap} ${refcells} \\
            ${group_id} numbat ${task.cpus} \\
            ${params.numbat_min_cells} ${params.numbat_t} ${params.numbat_max_iter} \\
            ${params.numbat_min_llr} ${params.numbat_min_ref_cells} \\
            ${alleleList.join(' ')}

        # Numbat's own diagnostic plots, when it produced any.
        find numbat/${group_id} -name "*.png" -exec cp {} figures/numbat/ \\; 2>/dev/null || true
        """

    stub:
        """
        mkdir -p numbat/${group_id} figures/numbat
        printf 'cell\\tbarcode\\tcell_id\\tsample\\tgroup\\tnumbat_call\\tp_cnv\\n' > numbat/${group_id}_numbat_calls.tsv
        printf 'TESTA_AACGTAACGTAACGTA-1\\tAACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\t${group_id}\\tmalignant\\t0.98\\n' >> numbat/${group_id}_numbat_calls.tsv
        printf 'sample\\tn_cells\\tmalignant\\tpct_malignant\\treference_basis\\nTESTA\\t1\\t1\\t100.0\\tref_hca\\n' > numbat/${group_id}_numbat_by_library.tsv
        printf 'stub: null-result note\\n' > numbat/${group_id}_numbat_nocalls.txt
        touch numbat/${group_id}/stub.txt
        """
}


// Numbat's own figures, drawn from a finished NUMBAT_RUN output directory.
//
// A SEPARATE process on purpose. Nextflow hashes a process's script text, so
// adding the plotting to NUMBAT_RUN would re-run the model fit for every
// sample -- 957 CPU-hours over 15 libraries on GBM-DFCI3, to produce pictures
// from numbers that are already on disk. run_numbat() writes its posteriors,
// consensus segments, phylogeny and mutation graph OUTSIDE its `if (plot)`
// branch, so `plot = FALSE` cost us the figures and none of the data; this
// reads that data back and draws them in minutes.
//
// Runs in the Numbat image because that is where the numbat R package and its
// plotting stack live (ggtree, ggraph, tidygraph, patchwork, scistreer).
// There is no quarto in that image, so this emits figure FILES rather than a
// notebook -- which is what STAGE_REPORT consumes anyway: it builds the CNV
// report from what the stage published.
process NUMBAT_REPORT {

    tag "Numbat figures: ${group_id}"
    label 'process_medium'

    input:
        tuple val(group_id), path(numbat_dir)

    output:
        path "figures/numbat/*",                  emit: figures,  optional: true
        path "data/*_numbat_figure_manifest.tsv", emit: manifest, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p figures/numbat data

        # The last two positional arguments are the rolling-expression cluster
        # count and the cell count above which the dense heatmaps are written
        # as PNG only. Neither is worth a pipeline parameter: the first is
        # Numbat's own default, and the second is a property of the PDF format
        # rather than of the analysis.
        numbat_figures.R \\
            ${numbat_dir} ${group_id} figures/numbat data \\
            ${params.numbat_genome} \\
            ${params.numbat_figures_ncut} \\
            ${params.numbat_figures_min_clone_cells} \\
            ${params.numbat_figures_min_llr} \\
            3 3000
        """

    stub:
        """
        mkdir -p figures/numbat data
        touch figures/numbat/${group_id}_phylo_heatmap.png
        touch figures/numbat/${group_id}_consensus_segments.png
        printf 'sample\\tfigure\\tstatus\\treason\\n${group_id}\\tphylo_heatmap\\tdrawn\\t\\n' \\
            > data/${group_id}_numbat_figure_manifest.tsv
        """
}
