process CELLRANGER_COUNT {

    tag "Running GEX on ${sample}"
    label 'process_high'


    input:
        tuple val(sample), path(reads)
        path(reference)

    output:
        tuple val(sample), path("${sample}/outs/*"), emit: outs

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
            cellranger \\
                count \\
                --id="${sample}" \\
                --fastqs=. \\
                --transcriptome="${reference.name}" \\
                --sample="${sample}" \\
                --localcores=${task.cpus} \\
                --localmem=${task.memory.toGiga()} \\
                --create-bam true


            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
                reference: ${reference.name}
            END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        """
            mkdir -p ${sample}/outs/filtered_feature_bc_matrix

            touch ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
            touch ${sample}/outs/filtered_feature_bc_matrix/features.tsv.gz
            touch ${sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
            touch ${sample}/outs/metrics_summary.csv

            # The QC stage keys off these two filenames specifically; the stub
            # must produce them or the downstream grouping yields nulls.
            touch ${sample}/outs/filtered_feature_bc_matrix.h5
            touch ${sample}/outs/possorted_genome_bam.bam

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
                reference: ${reference.name}
            END_VERSIONS
        """
}
