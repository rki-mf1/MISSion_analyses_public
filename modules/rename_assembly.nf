// rename Assembly
process rename_spades {

    label "rename_spades"

    tag "rename_spades_${sample}"

    executor 'local'

    input:
        tuple val(sample), path(contigs)
        tuple val(sample_scaff), path(scaffolds)
        val(method)

    output:
        tuple val(sample), path("${sample}_contigs.fasta"), emit: contigs
        tuple val(sample), path("${sample}_scaffolds.fasta"), emit: scaffolds

    script:
        """
            cp ${contigs} ${sample}_contigs.fasta
            cp ${scaffolds} ${sample}_scaffolds.fasta
        """

    stub:
        """
            touch ${sample}_contigs.fasta ${sample}_scaffolds.fasta
        """
}


process rename_skesa {

    label "rename_skesa"

    tag "rename_skesa_${sample}"

    executor 'local'

    publishDir path: "${params.analysesdir}/04_assembly/${method}/${sample}", failOnError: true, mode: "copy"

    input:
        tuple val(sample), path(contigs)
        tuple val(sample_graph), path(graph)
        val(method)

    output:
        tuple val(sample), path("${sample}_contigs.fasta"), emit: contigs
        tuple val(sample), path("${sample}_contigs.gfa"), emit: graph

    script:
        """
            cp -l ${contigs} ${sample}_contigs.fasta
            cp -l ${graph} ${sample}_contigs.gfa
        """

    stub:
        """
            touch ${sample}_contigs.fasta
            touch ${sample}_contigs.gfa
        """
}
