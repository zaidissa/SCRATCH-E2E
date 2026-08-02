process SAMPLESHEET_CHECK {
    tag "Samplesheet $samplesheet"
    label 'process_single'


    input:
        path samplesheet

    output:
        path 'samplesheet.valid.csv', emit: csv
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        cp ${samplesheet} samplesheet.valid.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """


    stub:
        // The stub must pass the real samplesheet through: everything downstream
        // is driven by splitCsv over this file, so emitting placeholder text
        // would starve the entire DAG and make -stub useless as a dry run.
        """
        cp ${samplesheet} samplesheet.valid.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            stub: true
        END_VERSIONS
        """
}