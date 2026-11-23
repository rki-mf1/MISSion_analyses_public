// this process may also be moved into a module, but is located here to show how "inline" processes are defined.
process trim {

    label "fastp"

    conda "${projectDir}/envs/fastp.yaml"

    // this line will publish all "fastq.gz" files to the output directory
    // publishDir "${params.analysesdir}/01_trimming/${sample}", pattern: "*.fastq.gz", mode: "copy"

    // if you want to publish all output that was defined in the "output:" section further down, one leaves out the "pattern:"
    // we're using this one, because we also want to publish the fastp trimming reports
    publishDir "${params.analysesdir}/01_trimming/${sample}", mode: "copy"

    // SLURM cluster options
    cpus 10
    memory "3 GB"
    time "1h"
    tag "trim_${sample}"

    input:
        tuple val(sample), path(read1), path(read2)

    // each output of fastp is tagged with an "emit" to directly reference it for downstream usage
    output:
        tuple val(sample), path("${sample}.R*.trimmed.fastq.gz"), emit: fastq
        tuple val(sample), path("${sample}.trim_report.json"), emit: json
        tuple val(sample), path("${sample}.trim_report.html"), emit: html


    // here we write what is going to be executed
    script:
        """
            fastp \
                --in1 ${read1} \
                --in2 ${read2} \
                --out1 ${sample}.R1.trimmed.fastq.gz \
                --out2 ${sample}.R2.trimmed.fastq.gz \
                --trim_poly_g \
                --json ${sample}.trim_report.json \
                --html ${sample}.trim_report.html \
                --thread ${task.cpus} \
                --report_title ${sample}
        """

    // we can start the pipeline in "--stub" mode. Then this stub is going to be executed that basically checks if the workflow processes fit together
    stub:
        """
            touch ${sample}.R{1,2}.trimmed.fastq.gz ${sample}.trim_report.json ${sample}.trim_report.html
        """

}
