process stats_compile {

    label "multiqc"

    conda "${projectDir}/envs/multiqc.yaml"

    publishDir "${params.analysesdir}/99_reporting/${method}_statistics", mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "25 GB"
    time "4h"

    tag "compile_mapping_statistics"

    input:
        path(reports, stageAs: 'temp??/*') 
        val(method)

    output:
        path("reporting_output/${method}_multiqc_report.html"),     emit: report 
        path("reporting_output/${method}_multiqc_report_data"),     emit: report_folder

    script:
        """
            multiqc \
                -o reporting_output/ \
                -n ${method}_multiqc_report.html \
                --title "MultiQC Report for ${method} data" \
                ${reports}

        """

    stub:
        """
            mkdir -p reporting_output
            touch reporting_output/${method}_multiqc_report.html
            mkdir -p reporting_output/${method}_multiqc_report_data
        """
}


