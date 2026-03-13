process MULTIQCMAPPER {
    tag 'multiqc_mapper'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "multiqc-mapper:0.1.0"

    input:
    path multiqc_data_dir
    path concepts_yaml

    output:
    path "unified.tsv", emit: unified
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc-mapper resolve \\
        $multiqc_data_dir \\
        $concepts_yaml \\
        --output unified.tsv \\
        --format tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc-mapper: \$(multiqc-mapper version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch unified.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc-mapper: \$(multiqc-mapper version)
    END_VERSIONS
    """
}
