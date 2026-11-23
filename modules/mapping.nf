process map_reads {

    label "map_reads"

    conda "${projectDir}/envs/bwa.yaml"


    // SLURM cluster options
    cpus 10
    memory "20 GB"
    time "3h"

    tag "map_${sample}"

    publishDir "${params.analysesdir}/03_mapping/${method}/${sample}", mode: "copy"

    input:
        tuple val(type), path(genome)
        tuple val(sample), path(reads)
        val(method)

    output:
        tuple val(sample), path("${method}/${sample}/${sample}.bam"), emit: bam
        tuple val(sample), path("${method}/${sample}/${sample}.sort.bam"), emit: bam_sorted
        tuple val(sample), path("${method}/${sample}/${sample}.sort.bam.bai"), emit: bai

    script:
        """

            INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
            mkdir -p ${method}/${sample}

            bwa-mem2 mem \\
                -t ${task.cpus} \\
                -R '@RG\\tID:${sample}\\tSM:${sample}' \\
                \$INDEX \\
                ${reads} | samtools view -@ ${task.cpus} -O BAM -o ${method}/${sample}/${sample}.bam

            samtools sort -O BAM -@ ${task.cpus} -o ${method}/${sample}/${sample}.sort.bam ${method}/${sample}/${sample}.bam

            samtools index ${method}/${sample}/${sample}.sort.bam

        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/${sample}.bam
            touch ${method}/${sample}/${sample}.sort.bam
            touch ${method}/${sample}/${sample}.sort.bam.bai
        """
}