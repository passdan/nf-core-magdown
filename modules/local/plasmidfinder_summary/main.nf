process PLASMIDFINDER_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04':
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(tsvs)

    output:
    tuple val(meta), path("*.plasmidfinder_summary.tsv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    out=${prefix}.plasmidfinder_summary.tsv

    # Build header dynamically from the first report so it tracks the
    # PlasmidFinder version, prefixed with sample + mag columns.
    first=\$(ls -1 *.tsv 2>/dev/null | grep -v "\${out}\$" | head -n 1)
    printf 'sample\\tmag\\t' > "\$out"
    head -n 1 "\$first" >> "\$out"

    # For each per-MAG report, drop the header and keep only data rows,
    # tagging each with its sample and mag id. MAGs with no hits add nothing.
    for f in \$(ls -1 *.tsv | grep -v "\${out}\$"); do
        mag=\$(basename "\$f" .tsv)
        tail -n +2 "\$f" | awk -v s="${meta.id}" -v m="\$mag" 'NF { print s"\\t"m"\\t"\$0 }' >> "\$out"
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plasmidfinder_summary.tsv
    """
}
