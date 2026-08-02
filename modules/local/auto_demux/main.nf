process AUTO_DEMUX {

    tag "Auto-demux on ${run_id}"
    label 'process_medium'


    input:
        // i2 is a string path (val) so that an empty value ('') is accepted
        // when the library was sequenced with single indexing — avoids staging
        // an empty path object.
        tuple val(run_id), path(r1), path(r2), path(i1), val(i2), val(modality)

    output:
        // demux_dir contains all per-sample FASTQs and is consumed by the
        // workflow's flatMap to create one channel element per sample.
        tuple path("demuxed/"), path("demuxed/detected_samples.csv"), emit: demux_result
        path("demuxed/demux_stats.csv"), emit: stats
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args      = task.ext.args      ?: ''
        def min_freq  = params.demux_min_freq   ?: 0.01
        def n_detect  = params.demux_n_detect   ?: 50000
        def mismatches= params.demux_mismatches ?: 1
        // Only add --i2 if the user provided a non-empty path
        def i2_arg    = (i2 && i2.trim()) ? "--i2 ${i2}" : ''
        """
        auto_demux.py \\
            --run-id    "${run_id}" \\
            --r1        ${r1} \\
            --r2        ${r2} \\
            --i1        ${i1} \\
            ${i2_arg} \\
            --modality  "${modality}" \\
            --n-detect  ${n_detect} \\
            --min-freq  ${min_freq} \\
            --mismatches ${mismatches} \\
            --outdir    demuxed \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
            auto_demux: 1.0.0
        END_VERSIONS
        """

    stub:
        """
        mkdir -p demuxed
        echo "sample,centroid,modality"                              > demuxed/detected_samples.csv
        echo "${run_id}_S01_ATCGATCG,ATCGATCG,${modality}"         >> demuxed/detected_samples.csv
        echo "${run_id}_S02_GCTAGCTA,GCTAGCTA,${modality}"         >> demuxed/detected_samples.csv

        # Stub FASTQs in cellranger naming format
        touch demuxed/${run_id}_S01_ATCGATCG_S1_L001_R1_001.fastq.gz
        touch demuxed/${run_id}_S01_ATCGATCG_S1_L001_R2_001.fastq.gz
        touch demuxed/${run_id}_S02_GCTAGCTA_S1_L001_R1_001.fastq.gz
        touch demuxed/${run_id}_S02_GCTAGCTA_S1_L001_R2_001.fastq.gz
        touch demuxed/Undetermined_S0_L001_R1_001.fastq.gz
        touch demuxed/Undetermined_S0_L001_R2_001.fastq.gz

        echo "sample,reads,fraction"                                 > demuxed/demux_stats.csv
        echo "${run_id}_S01_ATCGATCG,950000,0.9500"                >> demuxed/demux_stats.csv
        echo "${run_id}_S02_GCTAGCTA,2000,0.0020"                  >> demuxed/demux_stats.csv
        echo "Undetermined,48000,0.0480"                            >> demuxed/demux_stats.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //')
            auto_demux: 1.0.0
        END_VERSIONS
        """
}
