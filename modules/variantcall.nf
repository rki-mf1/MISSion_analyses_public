process variantcall_single {

    label "variantcall_freebayes"

    conda "${projectDir}/envs/freebayes.yaml"

    errorStrategy 'ignore'


    // SLURM cluster options
    cpus 1
    memory "20 GB"
    
    tag "variantcall_${sample}"

    publishDir "${params.analysesdir}/08_variantcall/${sample}", mode: "copy"

    input:
        tuple val(sample), path(bam)
        tuple val(sample_bai), path(bai)
        path(ref)
        val(method)
        
    output:
        tuple val(sample), path("${method}/${sample}/${sample}.g.vcf"), emit: variants
    
    script:
        """
            
            mkdir -p ${method}/${sample}
            
            freebayes \\
                -f ${ref} \\
                --haplotype-length -1 \\
                --gvcf \\
                ${bam} > ${method}/${sample}/${sample}.g.vcf
                
        """

    stub:
        """
            mkdir -p ${method}/${sample}
            touch ${method}/${sample}/${sample}.g.vcf
        """
}


process variantcall_multi {

    label "variantcall_freebayes"

    conda "${projectDir}/envs/freebayes.yaml"


    // SLURM cluster options
    time "90h"
    cpus 60
    memory "90 GB"
    
    tag "variantcall_all"

    publishDir "${params.analysesdir}/08_variantcall", mode: "copy"

    input:
        path(list_bam)
        path(list_bai)
        path(ref)
        val(method)
        
    output:
        path("${method}/all_samples.vcf.gz"), emit: variants
    
    script:
        """
            mkdir -p ${method}
            samtools faidx ${ref}
            
            #/usr/bin/time -v \\
            #freebayes \\
            #    -f ${ref} \\
            #    --haplotype-length -1 \\
            #    ${list_bam} | \\
            #    gzip > ${method}/all_samples.vcf.gz

            freebayes-parallel <(fasta_generate_regions.py ${ref}.fai 100000) \\
                ${task.cpus} \\
                -f ${ref} \\
                --haplotype-length -1 \\
                ${list_bam} | \\
                gzip > ${method}/all_samples.vcf.gz
                
                
        """

    stub:
        """
            mkdir -p ${method}
            touch ${method}/all_samples.vcf.gz
        """
}
