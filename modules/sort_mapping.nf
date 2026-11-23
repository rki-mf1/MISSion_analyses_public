process sort_mapping {

    label "sort_map"

    conda "${projectDir}/envs/samtools.yaml"


    // SLURM cluster options
    cpus 10
    memory "5 GB"
    time "3h"

    tag "map_sort_${sample}"

    input:
        tuple val(sample), path(bam)
        val(method)

    output:
        tuple val(sample), path("${method}/${sample}/${sample}.sort.bam"), emit: sort_bam
        
    script:
        """
            mkdir -p ${method}/${sample}
            samtools sort -@ ${task.cpu} \
            -o ${method}/${sample}/${sample}.sort.bam \
            ${bam}

        """

    stub:
        """
            mkdir -p ${method}/${sample}
            ${method}/${sample}/${sample}.sort.bam
        """
}