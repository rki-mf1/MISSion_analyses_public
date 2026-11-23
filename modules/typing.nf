process type_contigs_sistr {

    label "sistr"

    conda "${projectDir}/envs/sistr.yaml"
    errorStrategy 'ignore'


    publishDir "${params.analysesdir}/07_typing/${sample}", mode: "copy"

    // SLURM cluster options
    cpus 10
    memory "15 GB"
    time "4h"
    tag "typing_${sample}"

    input:
        tuple val(sample), path(contigs) 
        val(method)

    output:
        tuple val(sample), path("${sample}-allele-results.json"),    emit: sistr_alleles
        tuple val(sample), path("${sample}-novel-alleles.fasta"),    emit: sistr_novel_alleles
        tuple val(sample), path("${sample}-cgmlst-profile.csv"),     emit: sistr_cgmlst
        tuple val(sample), path("${sample}-sistr-output.tab"),       emit: sistr_output

    script:
        """
            mkdir -p tmp
            sistr \
                --qc \
                -vv \
                --alleles-output ${sample}-allele-results.json \
                --novel-alleles ${sample}-novel-alleles.fasta \
                --cgmlst-profiles ${sample}-cgmlst-profile.csv \
                --output-prediction ${sample}-sistr-output.tab \
                --output-format tab \
                --tmp-dir tmp \
                --threads ${task.cpus} \
                ${contigs}

        """

    stub:
        """
            touch ${sample}-allele-results.json \
                ${sample}-novel-alleles.fasta \
                ${sample}-cgmlst-profile.csv \
                ${sample}-sistr-output.tab
        """
}
