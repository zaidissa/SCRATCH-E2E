/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Numbat — allele-aware CNV calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ON by default (--run_numbat, true since a588a94). Unlike the other callers
    it needs the aligned BAM plus a phasing reference panel; a run with no BAM
    switches it off with a warning rather than failing.

    What it adds over inferCNV / SCEVAN / CopyKAT: those infer copy number from
    expression smoothed along the genome, so they are structurally blind to
    copy-neutral LOH — an event that changes allele ratios without changing
    dosage produces no expression signal. Numbat reads phased heterozygous SNPs
    from the BAM and sees it.

    Four processes, kept separate so the expensive pileup caches independently
    of the phasing and the model fit:

      NUMBAT_PREP     annotated object -> per-individual counts + an explicit
                      barcode -> cell_id map (the BAM speaks raw cellranger
                      barcodes; the merged object does not)
      NUMBAT_CELLSNP  BAM + barcodes  -> allele pileup            [per library]
      NUMBAT_PHASE    pileups -> joint genotype, Eagle2 phase,
                      phased allele counts                        [per individual]
      NUMBAT_RUN      counts + alleles -> per-cell calls          [per individual]

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
        path "numbat_input/*_refgroups.tsv",     emit: refgroups
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
            numbat_input ${helpers} \\
            ${params.numbat_cluster_col} \\
            ${params.numbat_ref_min_group_cells} \\
            '${params.numbat_extra_reference_clusters}' \\
            ${params.numbat_sample_type_col}
        """

    stub:
        """
        mkdir -p numbat_input
        touch numbat_input/PT1_counts.rds
        printf 'cell\\tcelltype\\tcluster\\tref_group\\n' > numbat_input/PT1_refcells.tsv
        printf 'ref_group\\tcelltype\\tcluster\\tn_cells\\tadmitted\\tbasis\\n' > numbat_input/PT1_refgroups.tsv
        echo "AACGTAACGTAACGTA-1" > numbat_input/TESTA_barcodes.tsv
        echo "AACGTAACGTAACGTC-1" > numbat_input/TESTB_barcodes.tsv
        printf 'cell_id\\tbarcode\\tcell\\tlibrary\\tgroup\\tin_reference_lineage\\tsample_type\\n' > numbat_input/PT1_cellmap.tsv
        printf 'TESTA_AACGTAACGTAACGTA-1\\tAACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\tPT1\\tFALSE\\tTP\\n' >> numbat_input/PT1_cellmap.tsv
        printf 'TESTB_AACGTAACGTAACGTC-1\\tAACGTAACGTAACGTC-1\\tTESTB_AACGTAACGTAACGTC-1\\tTESTB\\tPT1\\tTRUE\\tNB\\n' >> numbat_input/PT1_cellmap.tsv
        printf 'library\\tgroup\\nTESTA\\tPT1\\nTESTB\\tPT1\\n' > numbat_input/groups.tsv
        printf 'group\\tn_libraries\\tlibraries\\tn_cells\\tn_reference_cells\\tn_reference_usable\\tn_reference_columns\\treference_basis\\tsample_types\\nPT1\\t2\\tTESTA,TESTB\\t2\\t0\\t0\\t0\\tref_hca\\tNB,TP\\n' > numbat_input/grouping_report.tsv
        echo PT1 > numbat_input/samples.txt
        """
}


process NUMBAT_CELLSNP {

    tag "Numbat pileup: ${lib} (${group_id})"
    label 'process_high'

    // The heavy step: cellsnp-lite walks ONE library's BAM for allele counts at
    // the 1000G SNP positions. One task per LIBRARY, so an individual's
    // libraries pile up side by side.
    //
    // This is the first half of Numbat's pileup_and_phase.R, split off. That
    // script builds one cellsnp-lite command per sample and runs them one after
    // another (`sh run_pileup.sh`), so for one individual wall time was the SUM
    // over its libraries -- for GBM1-DFCI3, 17 libraries in one task. Nothing
    // in a pileup depends on another library: the joint step that makes
    // per-individual Numbat worth having is the genotyping and phasing, which
    // NUMBAT_PHASE still does over all of the individual's libraries at once.
    //
    // The command line is pileup_and_phase.R's 10x branch, flag for flag
    // (numbat-rbase:latest, bin/pileup_and_phase.R L163-174). The thread count
    // only changes how fast it runs, not what it writes, so the tables are the
    // ones the single-task version would have produced.

    input:
        tuple val(group_id), val(lib), path(bam), path(bai), path(barcodes)
        path  snpvcf

    output:
        // A directory named after the library. NUMBAT_PHASE puts it where
        // pileup_and_phase.R expects it: {outdir}/pileup/{sample}/.
        tuple val(group_id), val(lib), path("${lib}", type: 'dir'), emit: pileup

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        # A staged NO_FILE placeholder means "not supplied -- use the image copy".
        # The resolution and the checks are the same as in NUMBAT_PHASE.
        SNPVCF="${snpvcf}"
        case "\$SNPVCF" in *NO_FILE*) SNPVCF=\$(ls /data/genome1K.phase3.SNP_AF5e2.chr1toX.${params.numbat_genome}.vcf* 2>/dev/null | head -1) ;; esac
        [ -s "\$SNPVCF" ] || { echo "numbat SNP VCF '\$SNPVCF' is missing or empty." >&2; exit 1; }
        echo "[numbat] individual=${group_id} library=${lib} snpvcf=\$SNPVCF"

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
        if [ ! -s "${bam}.bai" ]; then
            echo "[numbat] no .bai staged for ${bam} -- indexing (slow)" >&2
            LD_LIBRARY_PATH=/lib/x86_64-linux-gnu samtools index -@ ${task.cpus} "${bam}"
        fi

        cellsnp-lite \\
            -s ${bam} \\
            -b ${barcodes} \\
            -O ${lib} \\
            -R "\$SNPVCF" \\
            -p ${task.cpus} \\
            --minMAF 0 \\
            --minCOUNT 2 \\
            --UMItag Auto \\
            --cellTAG CB

        # pileup_and_phase.R reads exactly these four files per sample. Check
        # them HERE, where a failure costs one library and retries alone, rather
        # than in NUMBAT_PHASE, where it would cost the individual. An empty VCF
        # is checked as well: upstream stops on "0 variants" too, but only after
        # every other library has been piled up.
        for f in cellSNP.base.vcf cellSNP.tag.AD.mtx cellSNP.tag.DP.mtx cellSNP.samples.tsv; do
            [ -s "${lib}/\$f" ] || { echo "cellsnp-lite wrote no ${lib}/\$f" >&2; exit 1; }
        done
        if ! grep -v -m1 '^#' "${lib}/cellSNP.base.vcf" >/dev/null; then
            echo "pileup VCF for ${lib} has 0 variants" >&2; exit 1
        fi

        # The OTH matrix is never read downstream. It is as large as AD and DP,
        # and every byte of a declared output is copied to and from S3 on Cirro.
        rm -f "${lib}/cellSNP.tag.OTH.mtx"
        """

    stub:
        """
        mkdir -p ${lib}
        for f in cellSNP.base.vcf cellSNP.tag.AD.mtx cellSNP.tag.DP.mtx cellSNP.samples.tsv; do
            echo stub > ${lib}/\$f
        done
        """
}


process NUMBAT_PHASE {

    tag "Numbat phasing: ${group_id}"
    label 'process_medium'

    // ONE TASK PER INDIVIDUAL: the second half of pileup_and_phase.R, over the
    // pileups NUMBAT_CELLSNP made one library at a time.
    //
    // This is the part that has to be joint. pileup_and_phase.R says so itself
    // ("All samples must belong to the same individual"): it genotypes the
    // libraries together (numbat:::genotype(label, samples, vcfs, ...)), so a
    // heterozygous site gets coverage from every library the individual
    // contributed -- the allelic power expression-only callers lack -- and
    // Eagle2 phases that one genotype. Each library's allele table is then cut
    // from the shared phase.
    //
    // It runs Numbat's own script, minus the pileup block
    // (bin/numbat_phase_only.sh says what is cut and refuses to run if the
    // upstream script no longer has that shape). No BAM is staged: the script
    // only checks that --bams and --barcodes EXIST before reaching the cut
    // block, and the real files are hundreds of GB on S3 for an individual
    // like GBM1-DFCI3. Empty placeholders satisfy those checks, and nothing
    // after the cut reads them.
    //
    // stageAs 'snp??/*': each library's pileup directory keeps its own name in
    // a numbered directory, so no library name can collide with anything the
    // script creates in the task directory. The subworkflow sorts the rows by
    // library, so the order -- and the task hash -- are stable across runs.

    input:
        tuple val(group_id), val(libs), path(pileups, stageAs: 'snp??/*')
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
        def libList = (libs    instanceof List ? libs    : [libs]).collect { it.toString() }
        def dirList = (pileups instanceof List ? pileups : [pileups]).collect { it.toString() }
        if (libList.size() != dirList.size())
            error "NUMBAT_PHASE ${group_id}: ${libList.size()} librar(ies) but " +
                  "${dirList.size()} pileup director(ies) -- these must match."
        """
        set -euo pipefail
        mkdir -p pileup/pileup placeholders

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
        numbat_phase_only.sh "\$PILEUP" phase_only.R

        # Each library's pileup goes where the script reads it,
        # {outdir}/pileup/{sample}/, and gets a placeholder BAM and barcode file
        # for the existence checks (see the header comment).
        libs=(${libList.join(' ')})
        dirs=(${dirList.join(' ')})
        bams=(); bcs=()
        for i in "\${!libs[@]}"; do
            lb="\${libs[\$i]}"; d="\${dirs[\$i]}"
            [ "\$(basename "\$d")" = "\$lb" ] || {
                echo "pileup directory \$d does not belong to library \$lb" >&2; exit 1; }
            ln -s "\$PWD/\$d" "pileup/pileup/\$lb"
            : > "placeholders/\$lb.bam"; : > "placeholders/\$lb.barcodes.tsv"
            bams+=("placeholders/\$lb.bam"); bcs+=("placeholders/\$lb.barcodes.tsv")
        done
        csv() { local IFS=,; echo "\$*"; }

        Rscript phase_only.R \\
            --label     ${group_id} \\
            --samples   "\$(csv "\${libs[@]}")" \\
            --bams      "\$(csv "\${bams[@]}")" \\
            --barcodes  "\$(csv "\${bcs[@]}")" \\
            --outdir    pileup \\
            --gmap      "\$GMAP" \\
            --snpvcf    "\$SNPVCF" \\
            --paneldir  "\$PANELDIR" \\
            --ncores    ${task.cpus}

        # pileup_and_phase.R writes {outdir}/{sample}_allele_counts.tsv.gz once
        # per library. Fail loudly if a library produced none: an allele table
        # that is quietly absent costs that library its allelic channel, which is
        # the only thing Numbat has over inferCNV, and nothing downstream can
        # tell the difference between "no table" and "no events".
        for lb in "\${libs[@]}"; do
            [ -s "pileup/\${lb}_allele_counts.tsv.gz" ] || {
                echo "phasing produced no allele table for library \$lb of ${group_id}." >&2
                exit 1; }
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
        // The same figures pooled by COMPARTMENT. "There should be almost zero
        // neoplastic cells in PBMC" is a claim about a compartment, and pooling
        // its libraries is also what stops a 358-cell library being read as the
        // peer of a 10,000-cell one. Optional: written only when the metadata
        // carried a sample-type column at all.
        path "numbat/${group_id}_numbat_by_sample_type.tsv", emit: by_sample_type, optional: true
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
        echo "[numbat] group=${group_id} libraries=${alleleList.size()} min_LLR=${params.numbat_min_llr} max_entropy=${params.numbat_max_entropy} min_cells=${params.numbat_min_cells} t=${params.numbat_t} max_iter=${params.numbat_max_iter}"

        numbat_run.R \\
            ${counts} ${cellmap} ${refcells} \\
            ${group_id} numbat ${task.cpus} \\
            ${params.numbat_min_cells} ${params.numbat_t} ${params.numbat_max_iter} \\
            ${params.numbat_min_llr} ${params.numbat_min_ref_cells} \\
            ${params.numbat_max_entropy} \\
            ${alleleList.join(' ')}

        # Numbat's own diagnostic plots, when it produced any.
        find numbat/${group_id} -name "*.png" -exec cp {} figures/numbat/ \\; 2>/dev/null || true
        """

    stub:
        """
        mkdir -p numbat/${group_id} figures/numbat
        printf 'cell\\tbarcode\\tcell_id\\tsample\\tgroup\\tnumbat_call\\tp_cnv\\n' > numbat/${group_id}_numbat_calls.tsv
        printf 'TESTA_AACGTAACGTAACGTA-1\\tAACGTAACGTAACGTA-1\\tTESTA_AACGTAACGTAACGTA-1\\tTESTA\\t${group_id}\\tmalignant\\t0.98\\n' >> numbat/${group_id}_numbat_calls.tsv
        printf 'sample\\tsample_type\\tn_cells\\tn_reference_lineage\\tmalignant\\tpct_malignant_raw\\tmalignant_gated\\tpct_malignant_gated\\tgate_overrode\\treference_basis\\tn_reference_columns\\nTESTA\\tTP\\t1\\t0\\t1\\t100.0\\t1\\t100.0\\t0\\tref_hca\\tNA\\n' > numbat/${group_id}_numbat_by_library.tsv
        printf 'sample_type\\tn_libraries\\tn_cells\\tn_reference_lineage\\tmalignant\\tpct_malignant_raw\\tmalignant_gated\\tpct_malignant_gated\\tgate_overrode\\treference_basis\\nTP\\t1\\t1\\t0\\t1\\t100.0\\t1\\t100.0\\t0\\tref_hca\\n' > numbat/${group_id}_numbat_by_sample_type.tsv
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
