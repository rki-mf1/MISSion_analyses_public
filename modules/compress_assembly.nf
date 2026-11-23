process compress_spades {

    label "compress_spades"

    scratch '$TMPDIR'

    conda "${projectDir}/envs/pigz.yaml"
    
    // SLURM cluster options
    memory "4 GB"

    tag "compress_spades_${sample}"

    publishDir path: "${params.analysesdir}/04_assembly/${method}/${sample}", pattern: "*.fasta.gz", failOnError: true, mode: "copy"

    input:
        tuple val(sample), path(contigs)
        tuple val(sample_scaff), path(scaffolds)
        val(method)

    output:
        tuple val(sample), path("${sample}_contigs.fasta.gz"), emit: contigs
        tuple val(sample), path("${sample}_scaffolds.fasta.gz"), emit: scaffolds

    script:
        """
            pigz -c \$(realpath ${contigs}) > ${sample}_contigs.fasta.gz
            pigz -c \$(realpath ${scaffolds}) > ${sample}_scaffolds.fasta.gz
        """

    stub:
        """
            touch ${sample}_contigs.fasta.gz ${sample}_scaffolds.fasta.gz
        """
}

process compress_skesa {

    label "compress_skesa"

    scratch '$TMPDIR'

    conda "${projectDir}/envs/pigz.yaml"
    
    // SLURM cluster options
    memory "4 GB"

    tag "compress_skesa_${sample}"

    publishDir path: "${params.analysesdir}/04_assembly/${method}/${sample}", pattern: "*.fasta.gz", failOnError: true, mode: "copy"

    input:
        tuple val(sample), path(contigs)
        val(method)

    output:
        tuple val(sample), path("${sample}_contigs.fasta.gz"), emit: contigs

    script:
        """
            pigz -c \$(realpath ${contigs}) > ${sample}_contigs.fasta.gz
        """

    stub:
        """
            touch ${sample}_contigs.fasta.gz
        """
}