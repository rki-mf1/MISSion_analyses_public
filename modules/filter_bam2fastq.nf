process bam2fastq {

    label "bam2fastq"

    conda "${projectDir}/envs/samtools.yaml"


    // SLURM cluster options
    cpus 10
    memory "20 GB"
    time "3h"

    tag "transform_${sample}"

    publishDir "${params.analysesdir}/04_filtered_reads/${sample}", mode: "copy"

    input:
        tuple val(sample), path(bam)
        val(method)

    output:
        tuple val(sample), path("${method}/${sample}/${sample}.R*.fastq"), emit: bam
        // tuple val(sample), path("${method}/${sample}/${sample}.bam.bai"), emit: bai

    script:
        """

            samtools fastq \
                -@ ${task.cpus} \
                -F 12 \
                -1 ${method}/${sample}.R1.fastq \
                -2 ${method}/${sample}.R2.fastq \
                ${bam}

            pigz -p ${task.cpus} ${method}/${sample}.R1.fastq
            pigz -p ${task.cpus} ${method}/${sample}.R2.fastq

        """

    stub:
        """
            mkdir -p ${method}
            touch ${method}/${sample}.R1.fastq.gz ${method}/${sample}.R2.fastq.gz
        """
}