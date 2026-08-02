process CELLRANGER_MULTI {

    tag "Running cellranger multi on ${sample_id}"
    label 'process_high'


    input:
        tuple val(sample_id), path(multi_config)

    output:
        tuple val(sample_id), path("${sample_id}/outs/per_sample_outs"), emit: per_sample_outs
        tuple val(sample_id), path("${sample_id}/outs/multi"),           emit: multi_outs
        path "versions.yml",                                              emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        cellranger multi \\
            --id="${sample_id}" \\
            --csv="${multi_config}" \\
            --localcores=${task.cpus} \\
            --localmem=${task.memory.toGiga()}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(echo \$(cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
        END_VERSIONS
        """

    stub:
        """
        mkdir -p ${sample_id}/outs/per_sample_outs/SAMPLE1/count/sample_filtered_feature_bc_matrix
        touch ${sample_id}/outs/per_sample_outs/SAMPLE1/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
        touch ${sample_id}/outs/per_sample_outs/SAMPLE1/count/sample_filtered_feature_bc_matrix/features.tsv.gz
        touch ${sample_id}/outs/per_sample_outs/SAMPLE1/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz
        touch ${sample_id}/outs/per_sample_outs/SAMPLE1/metrics_summary.csv
        mkdir -p ${sample_id}/outs/multi/multiplexing_analysis

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(echo \$(cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
        END_VERSIONS
        """
}
