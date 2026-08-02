process CELLRANGER_VDJ {

    tag "Running VDJ on ${sample}"
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
            cellranger_renaming.py "${sample}" .

            cellranger \\
                vdj \\
                --id="${sample}" \\
                --fastqs=. \\
                --reference="${reference.name}" \\
                --chain=TR \\
                --localcores=${task.cpus} \\
                --localmem=${task.memory.toGiga()}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
                reference: ${reference.name}
            END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        """
            cellranger_renaming.py "${sample}" .

            mkdir -p ${sample}/outs/
            touch ${sample}/outs/filtered_contig.fasta
            touch ${sample}/outs/filtered_contig_annotations.csv
            touch ${sample}/outs/filtered_contig.fastq
            touch ${sample}/outs/metrics_summary.csv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
                reference: ${reference.name}
            END_VERSIONS
        """
}
