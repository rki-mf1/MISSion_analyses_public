// index reference genome
process index_genome {

    label "index_genome"

    conda "${projectDir}/envs/bwa.yaml"
    //errorStrategy = 'ignore' // large assemblies can kill all the pipeline on indexing


    // SLURM cluster options
    cpus 1
    memory "50 GB"
    time "3h"

    tag "index_genome"

    input:
        tuple val(sample), path(reference)

    output:
        tuple val(sample), path("genome"), emit: index
        
    script:
        """
            mkdir -p genome
            bwa-mem2 index -p genome/genome ${reference}

        """

    stub:
        """
            mkdir -p genome
            touch genome/genome.0123
            touch genome/genome.amb
            touch genome/genome.ann
            touch genome/genome.bwt.2bit.64
            touch genome/genome.pac

        """
}


// index reference genome
process index_genome_big {

    label "index_genome"

    conda "${projectDir}/envs/bwa.yaml"


    // SLURM cluster options
    cpus 1
    memory "250 GB"
    time "20h"

    tag "index_genome"

    input:
        tuple val(sample), path(reference)

    output:
        tuple val(sample), path("genome"), emit: index
        
    script:
        """
            mkdir -p genome
            bwa-mem2 index -p genome/genome ${reference}

        """

    stub:
        """
            mkdir -p genome
            touch genome/genome.0123
            touch genome/genome.amb
            touch genome/genome.ann
            touch genome/genome.bwt.2bit.64
            touch genome/genome.pac

        """
}