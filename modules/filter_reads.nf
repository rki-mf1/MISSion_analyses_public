process filter_reads {

    label "krakentools"

    conda "${projectDir}/envs/krakentools.yaml"
    errorStrategy 'ignore'

    // SLURM cluster options
    cpus 1
    memory "15 GB"
    time "3h"
    // clusterOptions "--gres=local:50"

    tag "filter_${sample}"

    input:
        tuple val(sample), path(reads)
        tuple val(sample_out), path(kraken_output)
        tuple val(sample_report), path(kraken_report)
        val(taxid)

    output:
        tuple val(sample), file("${sample}.classified.R*.fastq"),     emit: fastq
    
    script:
        """
            extract_kraken_reads.py \
                -t ${taxid}\
                -k ${kraken_output} \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${sample}.classified.R1.fastq \
                -o2 ${sample}.classified.R2.fastq \
                --fastq-output \
                --report ${kraken_report} \
                --include-children

        """
    stub:
        """
            touch ${sample}.classified.R1.fastq ${sample}.classified.R2.fastq
        """
}

process filter_contigs {

    label "krakentools"

    conda "${projectDir}/envs/krakentools.yaml"
    
    errorStrategy 'ignore'

    // SLURM cluster options
    cpus 1
    memory "15 GB"
    time "3h"

    publishDir "${params.analysesdir}/04_assembly/${method}/${sample}", pattern: "*.fasta", mode: "copy"

    tag "filter_${sample}"

    input:
        tuple val(sample), path(contigs)
        tuple val(sample_out), path(kraken_output)
        tuple val(sample_report), path(kraken_report)
        val(taxid)
        val(method)

    output:
        tuple val(sample), file("${sample}.classified.fasta"),     emit: fastq
    
    script:
        """
            extract_kraken_reads.py \
                -t ${taxid}\
                -k ${kraken_output} \
                -s ${contigs} \
                -o ${sample}.classified.fasta \
                --report ${kraken_report} \
                --include-children

        """
    stub:
        """
            touch ${sample}.classified.fasta 
        """
}


process filter_out_reads {

    label "krakentools"

    conda "${projectDir}/envs/krakentools.yaml"

    // SLURM cluster options
    cpus 1
    memory "25 GB"
    time "5h"

    tag "filter_${sample}"

    input:
        tuple val(sample), path(reads)
        tuple val(sample_out), path(kraken_output)
        tuple val(sample_report), path(kraken_report)
        val(taxid)

    output:
        tuple val(sample), file("${sample}.classified.R*.fastq"),     emit: fastq
    
    script:
        """
            extract_kraken_reads.py \
                -t ${taxid}\
                -k ${kraken_output} \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${sample}.classified.R1.fastq \
                -o2 ${sample}.classified.R2.fastq \
                --fastq-output \
                --report ${kraken_report} \
                --include-children \
                --exclude

        """
    stub:
        """
            touch ${sample}.classified.R1.fastq ${sample}.classified.R2.fastq
        """
}


process filter_bam {

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
        tuple val(sample), path("${method}/${sample}.R*.fastq.gz"), emit: bam
        // tuple val(sample), path("${method}/${sample}/${sample}.bam.bai"), emit: bai

    script:
        """
            mkdir -p ${method}

            samtools fastq \
                -@ ${task.cpus} \
                -f 12 \
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