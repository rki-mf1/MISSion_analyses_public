process variantcall_single_statistics {

    label "variantcall_statistics"

    conda "${projectDir}/envs/bcftools.yaml"


    // SLURM cluster options
    cpus 1
    memory "5 GB"
    
    tag "variantstats_${sample}"

    publishDir "${params.analysesdir}/99_reporting/vcf_statistics/${sample}", mode: "copy"

    input:
        tuple val(sample), path(vcf)
                
    output:
        tuple val(sample), path("${sample}.csv"), emit: variantcounts
    
    script:
        """
            # get "heterozygous" aka mixed variants
            HET=\$(bcftools \\
                view \\
                -H \\
                -c1 \\
                --genotype 'het' ${vcf} | wc -l )

            HOM=\$(bcftools \\
                view \\
                -H \\
                -c1 \\
                --genotype 'hom' ${vcf} | wc -l )
            
            FRAC=\$(awk \\
                    -v het="\$HET" \\
                    -v hom="\$HOM" \\
                    'BEGIN {if (het+hom == 0){ print 0 } else {print het/(het+hom)}}')

            echo "SAMPLE,HET,HOM,HET_FRAC"      > ${sample}.csv
            echo "${sample},\$HET,\$HOM,\$FRAC" >> ${sample}.csv
        """

    stub:
        """
            touch ${sample}.csv
        """
}

process compile_variant_counts {

    label "compile_variant_counts"


    publishDir path: "${params.analysesdir}/99_reporting", failOnError: true, mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "1 GB"
    time "1h"

    tag "variant_count_compilation"

    input:
        path(count_files)

    output:
        path("variant_count_compiled.csv"), emit: compiled_variant_counts
        
    script:
        """
            echo "sample,het,hom,het_fraction" > variant_count_compiled.csv

            tail -q -n1 ${count_files} >> variant_count_compiled.csv
        """

    stub:
        """
            touch variant_count_compiled.csv
        """
}

