process HELPER_SCEASY_CONVERTER {

    tag "Converting ${seurat_object.baseName} to AnnData"
    label 'process_medium'

    

    input:
        path(seurat_object)

    output:
        path("${seurat_object.baseName}.h5ad") ,  emit: project_rds

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        sceasy_converter.R -f ${seurat_object}
        """
    stub:
        """
        touch ${seurat_object.baseName}.h5ad
        """

}
