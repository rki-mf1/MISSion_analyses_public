process qc_assembly {

    label "qc_assembly"

    conda "${projectDir}/envs/quast.yaml"

    // publishDir path: "${params.analysesdir}/05_qc_assembly", pattern: "${method}/${sample}/report.*", failOnError: true, mode: "copy"
    // publishDir path: "${params.analysesdir}/05_qc_assembly", pattern: "${method}/${sample}/icarus.*.html", failOnError: true, mode: "copy"
    publishDir path: "${params.analysesdir}/05_qc_assembly", failOnError: true, mode: "copy"

    // SLURM cluster options
    cpus 10
    memory "30 GB"
    time "1h"
    // clusterOptions "--job-name=assembly_${sample}"
    errorStrategy 'ignore'

    tag "assemblyqc_${sample}"

    input:
        tuple val(sample), path(contigs)
        val(fasta)
        val(gff)
        val(method)

    output:
        tuple val(sample), path("${sample}/${method}/report.html"), emit: report_html
        tuple val(sample), path("${sample}/${method}/report.tsv"), emit: report_tsv // this is the file multiQC needs. has to be named report.tsv
        tuple val(sample), path("${sample}/${method}/icarus.html"), emit: icarus_html

    script:
        """
            # version 4.6.0
            quast.py --output-dir ${sample} \
                --labels "${sample}_${method}" \
                -R ${fasta} \
                --features ${gff} \
                -t ${task.cpus} \
                ${contigs} 

            mkdir -p ${sample}/${method}
            ln -sr ${sample}/report.html ${sample}/${method}/report.html
            ln -sr ${sample}/report.tsv ${sample}/${method}/report.tsv
            ln -sr ${sample}/icarus.html ${sample}/${method}/icarus.html
        """

    stub:
        """
            mkdir -p ${sample}/${method}
            touch ${sample}/${method}/report.html
            touch ${sample}/${method}/report.tsv
            touch ${sample}/${method}/icarus.html
        """
}
