process map_stats {

    label "map_statistics"

    conda "${projectDir}/envs/samtools.yaml"


    // SLURM cluster options
    cpus 1
    memory "1 GB"
    time "3h"

    tag "map_stats_${sample}"

    publishDir "${params.analysesdir}/03_mapping/${sample}", mode: "copy"

    input:
        tuple val(sample), path(bam)
        val(method)

    output:
        tuple val(sample), path("${method}/${sample}/${sample}.flagstat"), emit: flagstat
        

    script:
        """
            mkdir -p ${method}/${sample}
            samtools flagstat ${bam} > ${method}/${sample}/${sample}.flagstat

        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/${sample}.flagstat
        """
}