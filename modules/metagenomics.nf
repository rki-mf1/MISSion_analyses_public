process metagenomics_setup {

    label "metaphlan_setup"

    conda "${projectDir}/envs/metaphlan.yaml"

    // SLURM cluster options
    cpus 11
    // memory "750 GB"
    memory { workflow.stubRun ? 1.GB : 4.GB}
    time "5"
    // clusterOptions "--job-name=classify_${sample}"
    tag "setup_metaphlan"

    output:
        path("setup.done.txt")

    script:
        """
            metaphlan \\
                --install \\
                --bowtie2db ${params.analysesdir}/09_metagenomes/databases && \\
            touch setup.done.txt
        """

    stub:
        """
            touch setup.done.txt
        """
}


process metagenomics {

    label "metaphlan"

    conda "${projectDir}/envs/metaphlan.yaml"

    publishDir "${params.analysesdir}/09_metagenomes/${method}/${sample}", mode: "copy"

    errorStrategy 'ignore'

    scratch true

    // SLURM cluster options
    cpus 10
    memory { workflow.stubRun ? 1.GB : 150.GB}
    time "24h"

    clusterOptions "--gres=local:250"
    
    tag "metagenome_${sample}"

    input:
        tuple val(sample), path(reads)
        val method
        val(bowtie2db)

    output:
        tuple val(sample), path("${sample}.txt"),     emit: profile
        tuple val(sample), path("${sample}.bowtie.bz2"), emit: profile_bowtie2

    script:
        """
            metaphlan \\
                ${reads[0]},${reads[1]} \\
                --input_type fastq \\
                --tmp_dir \$TMPDIR \\
                --nproc ${task.cpus} \\
                --bowtie2db ${bowtie2db} \\
                -o ${sample}.txt \\
                -x ${bowtie2db} \\
                --bowtie2out ${sample}.bowtie.bz2
                
        """

    stub:
        """
            touch ${sample}.txt ${sample}.bowtie.bz2
        """
}