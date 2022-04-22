process R {
    tag "$meta.sample_id"
    label 'process_medium'

    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*')

    script:
    def args = task.ext.args  ?: ''

    """
    Rscript ${params.script} \\
        --cores ${task.cpus} \\
        --runtype nextflow \\
        $args
        
    rm -r input
    """
}