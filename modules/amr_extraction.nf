process amr_extract {

    label "amr_extract"


    conda "${projectDir}/envs/seqtk.yaml"
    
    publishDir "${params.analysesdir}/06_abricate/${sample}", pattern: "*.fasta", mode: "copy"

    // SLURM cluster options
    cpus 1
    memory "5 GB"
    time "1h"

    tag "amr_extract_${sample}"

    input:
        tuple val(sample), path(amrfile), path(assembly)
        val(abricatedb)
        val(method)

    output:
        tuple val(sample), path("${sample}.${abricatedb}.${method}.fasta"), emit: amrgenes

    script:
        """
            # extract AMR gene positions from abricate output
            cut -f 2-4,6 ${amrfile} | tail -n+2 >> amr.bed

            # extract AMR gene sequences from assembly
            # watch out start position is offset by 1
            seqtk subseq ${assembly} amr.bed >> ${sample}.${abricatedb}.${method}.fasta

            # add gene name to fasta entry
            while read line
            do 
                contig=\$(echo \$line | cut -f1 -d' '); 
                start=\$(echo \$line | cut -f2 -d' '); 
                end=\$(echo \$line | cut -f3 -d ' '); 
                gene=\$(echo \$line | cut -f4 -d' ');

                start=\$(echo "\$start + 1" | bc )

                 sed -i -e "s@\${contig}:\${start}-\${end}@\${contig}:\${start}-\${end} \${gene}@" ${sample}.${abricatedb}.${method}.fasta; 
            done < amr.bed
        """

    stub:
        """
            touch ${sample}.${abricatedb}.${method}.fasta
        """
}
