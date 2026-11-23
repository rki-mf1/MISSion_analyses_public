process abricate {

    label "abricate"

    conda "${projectDir}/envs/abricate.yaml"
    // conda "/scratch/Projekte/MF1_Mission_kemri/tools/abricate_env"
    
    publishDir "${params.analysesdir}/06_abricate/${sample}", pattern: "*.tsv", mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "5 GB"
    time "4h"

    tag "abricate_${sample}"

    input:
        tuple val(sample), path(assembly)
        val(abricatedb)
        val(method)

    output:
        tuple val(sample), path("${sample}.${abricatedb}.${method}.tsv"), emit: amrlist

    script:
        """
            abricate \
                --db ${abricatedb} \
                ${assembly} \
                1> ${sample}.${abricatedb}.${method}.tsv
                2> ${sample}.${abricatedb}.${method}.log
        """

    stub:
        """
            touch ${sample}.${abricatedb}.${method}.tsv
        """
}