process classify {

    label "kraken"

    conda "${projectDir}/envs/kraken2.yaml"

    publishDir "${params.analysesdir}/02_classification/${method}/${sample}", pattern: "*.txt", mode: "copy"

    errorStrategy 'ignore'

    // SLURM cluster options
    cpus 10
    // memory "750 GB"
    memory { workflow.stubRun ? 1.GB : 750.GB}
    time "4h"
    // clusterOptions "--job-name=classify_${sample}"
    tag "classify_${sample}"

    input:
        // tuple val(sample), path(read1), path(read2) // <- if trimmed input from csv is used (trimming run with QCurchin)
        tuple val(sample), path(reads) // <- if fastp trimming process is used
        
        val krakendb
        val method

    output:
        tuple val(sample), path("${sample}.classified.R*.fastq"),     emit: fastq
        tuple val(sample), path("${sample}.kraken.out.txt"),          emit: kraken_output
        tuple val(sample), path("${sample}.kraken.report.txt"),       emit: kraken_report

    script:
        """
            kraken2 \
                --threads ${task.cpus} \
                --db ${krakendb} \
                --paired \
                --classified-out ${sample}.classified.R#.fastq \
                --output ${sample}.kraken.out.txt \
                --report ${sample}.kraken.report.txt \
                ${reads[0]} ${reads[1]}
        """

    stub:
        """
            touch ${sample}.classified.R_{1,2}.fastq ${sample}.kraken.out.txt ${sample}.kraken.report.txt
        """
}


process classify_contigs {

    label "kraken"

    conda "${projectDir}/envs/kraken2.yaml"

    publishDir "${params.analysesdir}/05_qc_assembly/${sample}", pattern: "*.txt", mode: "copy"
    publishDir "${params.analysesdir}/04_assembly/${method}/${sample}", pattern: "*.fasta", mode: "copy"

    // SLURM cluster options
    cpus 10
    memory { workflow.stubRun ? 1.GB : 750.GB}
    time "24h"

    tag "classify_${sample}"

    input:
        tuple val(sample), path(contigs) 
               
        val krakendb_contigs
        val method

    output:
        tuple val(sample), path("${sample}.classified_contigs.fasta"), emit: fasta
        tuple val(sample), path("${sample}.kraken.out.txt"),           emit: kraken_output
        tuple val(sample), path("${sample}.kraken.report.txt"),        emit: kraken_report

    script:
        """
            # gunzip -dc ${contigs} > ${sample}.fasta
            kraken2 \\
                --threads ${task.cpus} \\
                --db ${krakendb_contigs} \\
                --classified-out ${sample}.classified_contigs.fasta \\
                --output ${sample}.kraken.out.txt \\
                --report ${sample}.kraken.report.txt \\
                ${contigs}
        """

    stub:
        """
            touch ${sample}.classified_contigs.fasta ${sample}.kraken.out.txt ${sample}.kraken.report.txt
        """
}

process abundance {

    label "bracken"

    conda "${projectDir}/envs/bracken.yaml"

    publishDir "${params.analysesdir}/02_classification/${sample}", mode: "copy"

    errorStrategy 'ignore'

    // SLURM cluster options
    cpus 1
    memory { workflow.stubRun ? 1.GB : 750.GB}
    time "4h"
    
    tag "abundance_${sample}"

    input:
        tuple val(sample), path(report)
        val krakendb
        val level

    output:
        tuple val(sample), path("${sample}.${level}.bracken"), emit: bracken
        tuple val(sample), path("${sample}.${level}.bracken.report"), emit: bracken_report

    script:
        """
            # use defaults for readlength, level and threshold
            bracken \
                -d ${krakendb} \
                -i ${report} \
                -o ${sample}.${level}.bracken \
                -w ${sample}.${level}.bracken.report \
                -r 100 \
                -l ${level} \
                -t 0
        """

    stub:
        """
            touch ${sample}.${level}.bracken \
                  ${sample}.${level}.brakcen.report
        """
}

