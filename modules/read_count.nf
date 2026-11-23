process count_reads {

    label "count_reads"

    publishDir path: "${params.analysesdir}/${stage}/${sample}/", pattern: "*.csv", mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "15 GB"
    time "1h"

    tag "readcount_${sample}"

    input:
        tuple val(sample), path(reads)
        val(stage)
        
    output:
        tuple val(sample), path("${sample}_read_count.csv"), emit: read_count
        
    script:
        """
            zcat ${reads[0]} | awk 'END {printf "%s", "${sample},"; printf "%.0f", NR/4; print ""}' > ${sample}_read_count.csv

        """

    stub:
        """
            mkdir -p ${sample}
            touch ${sample}_read_count.csv
        """
}


process compile_read_counts {

    label "compile_read_counts"


    publishDir path: "${params.analysesdir}/99_reporting", failOnError: true, mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "1 GB"
    time "1h"

    tag "read_count_compilation"

    input:
        path(count_files)
        val(method)

    output:
        path("${method}_read_count_compiled.csv"), emit: compiled_read_counts
        
    script:
        """
            echo "sample,read_count" > ${method}_read_count_compiled.csv

            cat ${count_files} >> ${method}_read_count_compiled.csv
        """

    stub:
        """
            touch ${method}_read_count_compiled.csv
        """
}
