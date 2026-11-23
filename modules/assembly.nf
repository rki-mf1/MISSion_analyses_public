process assemble_spades {

    label "assemble_spades"

    conda "${projectDir}/envs/spades.yaml"
    errorStrategy 'ignore'


    // SLURM cluster options
    cpus 30
    // Spades fails, if there are too few reads. hence errorStrategy = 'ignore'
    // but
    // error code that Spades emits for out of memory - this should be compensated by increasing memory allocation
    // memory = { task.attempt < 4 ? 20.GB * task.attempt : 250.GB } ; maxRetries = 4 ; errorStrategy = { task.exitStatus == 247 ? 'retry' : 'ignore' } 
    memory "60 GB"
    
    tag "assemble_${sample}"

    input:
        tuple val(sample), path(reads)
        val(method)
        val(spades_mode)

    output:
        tuple val(sample), path("${method}/${sample}/contigs.fasta"), emit: contigs
        tuple val(sample), path("${method}/${sample}/scaffolds.fasta"), emit: scaffolds

    script:
        """
            spades.py \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                -t ${task.cpus} \
                ${spades_mode} \
                -o ${method}/${sample} &> ${sample}.spades.log
                
        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/contigs.fasta ${method}/${sample}/scaffolds.fasta
        """
}


process assemble_shovill {

    label "assemble_shovill"

    conda "${projectDir}/envs/shovill.yaml"
    errorStrategy 'ignore'


    // SLURM cluster options
    cpus 20
    memory "60 GB"
    
    tag "assemble_${sample}"

    input:
        tuple val(sample), path(reads)
        val(method)

    output:
        tuple val(sample), path("${method}/${sample}/contigs.fa"),  emit: contigs
        tuple val(sample), path("${method}/${sample}/spades.gfa"), emit: graph

    script:
        """
            shovill \
                --R1 ${reads[0]} \
                --R2 ${reads[1]} \
                --gsize 5M \
                --cpus ${task.cpus} \
                --ram ${task.memory} \
                --outdir ${method}/${sample} &> ${sample}.shovill.log
                
        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/contigs.fa
            touch ${method}/${sample}/spades.gfa
        """
}

process assemble_spades_meta {

    label "assemble_spades_meta"

    conda "${projectDir}/envs/spades.yaml"

    errorStrategy 'ignore'
    cpus 30
    // Spades fails, if there are too few reads. hence errorStrategy = 'ignore'
    // but
    // error code that Spades emits for out of memory - this should be compensated by increasing memory allocation
    // memory = { task.attempt < 4 ? 20.GB * task.attempt : 250.GB } ; maxRetries = 4 ; errorStrategy = { task.exitStatus == 247 ? 'retry' : 'ignore' } 
    memory { workflow.stubRun ? 1.GB : 500.GB}
    
    tag "assemble_meta_${sample}"

    publishDir "${params.analysesdir}/04_assembly/${method}/${sample}", mode: "copy"


    input:
        tuple val(sample), path(reads)
        val(method)
        val(spades_mode)

    output:
        tuple val(sample), path("${method}/${sample}/contigs.fasta"), emit: contigs
        tuple val(sample), path("${method}/${sample}/scaffolds.fasta"), emit: scaffolds //SRR5330172.contigs.gfa
        tuple val(sample), path("${method}/${sample}/contigs.gfa"), emit: graph

    script:
        """
            spades.py \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                -t ${task.cpus} \
                ${spades_mode} \
                -o ${method}/${sample} &> ${sample}.spades.log
                
        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/contigs.fasta ${method}/${sample}/scaffolds.fasta ${method}/${sample}/contigs.gfa
        """
}

