process CELLRANGER_MKFASTQ {

    tag "Running mkfastq on ${sample}"
    label 'process_high'



    input:
        tuple val(sample), val(run_dir), val(lane), val(index), val(modality)

    output:
        tuple val(sample), path("${sample}/outs/fastq_path"), val(modality), emit: fastq_dir
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        echo "Lane,Sample,Index" > ${sample}_simple.csv
        echo "${lane},${sample},${index}" >> ${sample}_simple.csv

        cellranger mkfastq \\
            --run="${run_dir}" \\
            --id="${sample}" \\
            --csv="${sample}_simple.csv" \\
            --localcores=${task.cpus} \\
            --localmem=${task.memory.toGiga()} \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
        END_VERSIONS
        """

    stub:
        """
        mkdir -p ${sample}/outs/fastq_path/STUBFLOWCELL/${sample}
        touch ${sample}/outs/fastq_path/STUBFLOWCELL/${sample}/${sample}_S1_L001_R1_001.fastq.gz
        touch ${sample}/outs/fastq_path/STUBFLOWCELL/${sample}/${sample}_S1_L001_R2_001.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
        END_VERSIONS
        """
}
