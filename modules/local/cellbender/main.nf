process CELLBENDER {

    tag "Running Cellbender on ${sample_id}"
    label 'process_high'

    

    input:
        tuple val(sample_id), path(csv_metrics), path(matrices)

    output:
        tuple val(sample_id), path(csv_metrics), path("cellbender_${sample_id}_matrix.h5")

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ? task.ext.args : ''
        def cellbender_output = "cellbender_${sample_id}_matrix.h5"
        """
        cellbender remove-background \
            --input ${matrices} \
            --output ${cellbender_output} \
            ${args}
        """
    stub:
        def args = task.ext.args ? task.ext.args : ''
        def cellbender_output = "cellbender_${sample_id}_matrix.h5"
        """
        touch ${cellbender_output}
        """

}