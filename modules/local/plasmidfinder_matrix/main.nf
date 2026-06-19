process PLASMIDFINDER_MATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2':
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(summaries)

    output:
    tuple val(meta), path("*.plasmidfinder_matrix.tsv"), emit: matrix
    tuple val(meta), path("*.plasmidfinder_long.tsv")  , emit: long

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plasmidfinder_matrix.py --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plasmidfinder_matrix.tsv
    touch ${prefix}.plasmidfinder_long.tsv
    """
}
