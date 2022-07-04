process MERGE_LOOM {
    tag "$meta.sample_id"
    label 'process_low'

    container "alexthiery/10x-npb-scvelo:base-2.0.1"

    input:
    tuple val(meta), path(loom_dir)

    output:
    tuple val(meta), path("*.loom"), emit: loom

    script:
    def prefix = task.ext.prefix ?: "merged"

    """
    $moduleDir/bin/merge_loom.py --input ${loom_dir} --output ${prefix}.loom
    """
}