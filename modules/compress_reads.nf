process compress_reads {

    label "compress"

    scratch '$TMPDIR'

    conda "${projectDir}/envs/pigz.yaml"

    publishDir path: "${params.analysesdir}/02_read_filter/${sample}", pattern: "*.fastq.gz", failOnError: true, mode: "copy"

    // SLURM cluster options
    cpus 10
    memory "2 GB"
    time "3h"

    tag "compress_${sample}"

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("${sample}.classified.R*.fastq.gz"), emit: fastq

    script:
        """
            pigz -p ${task.cpus} -c \$(realpath ${sample}.classified.R1.fastq) > ${sample}.classified.R1.fastq.gz
            pigz -p ${task.cpus} -c \$(realpath ${sample}.classified.R2.fastq) > ${sample}.classified.R2.fastq.gz
        """

    stub:
        """
            touch ${sample}.classified.R1.fastq.gz ${sample}.classified.R2.fastq.gz
        """

}
