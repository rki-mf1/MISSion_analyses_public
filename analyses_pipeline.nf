/************************** 
* CALL
**************************/

// (nextflow_env) 
// $ nextflow run analyses_pipeline.nf \
// --config samplesheet.csv \
// --reference_gff genomic.gff \
// --reference_fasta genomic.fna \
// --analysesdir results \
// -profile conda \
// -with-report logs/test.$(date +%T).report.html \
// #optional --krakendb references/kraken_db/ \ # pre-made virus DB for Kraken
// #optional --taxid 10239 \ # viruses
// #optional --krakendb_contigs \ # a krakenDB that is used to run taxonomic classification of reads after species filter (as QC) and classify contigs
// #optional --typing \ # run SISTR Salmonella strain typing
// #optional -resume
// #optional -stub



/************************** 
* DEFINE VARIABLES
**************************/


inputCSV      = params.config
reference_gff   = file(params.reference_gff)
reference_fasta = file(params.reference_fasta)
analysesdir   = params.analysesdir

if (params.krakendb && params.taxid) {
    krakendb      = params.krakendb // this is only necessary, if we want to do species filtering on environmental samples
    taxid         = params.taxid // this is only necessary, if we want to do species filtering on environmental samples, e.g. 10239 for viruses
    
    // check if another database is given for contig classification. If not, use the one for read filtering.
    krakendb_contigs = params.krakendb_contigs ? params.krakendb_contigs : params.krakendb 
    }
else {
    print("no Kraken db or taxid given in parameters. Switching to non-filtered assembly")
}





/************************** 
* CHECK VARIABLES
**************************/

// better safe than sorry, check parameter content
if(!params.config) {
    error("no sample sheet was defined.")
}
if(!params.reference_gff) {
    error("no GFF file was defined.")
}
if(!params.reference_fasta) {
    error("no FASTA file was defined.")
}

/************************** 
* DEFINE INCLUDES
**************************/
// necessary, because NextFlow does not allow to use processes twice, if not stored as a module

include { index_genome as index_reference } from "./modules/index_ref.nf"
include { index_genome as index_assembly}   from "./modules/index_ref.nf"
include { index_genome as index_assembly_meta}   from "./modules/index_ref.nf"
include { index_genome as index_mlst}   from "./modules/index_ref.nf"

// include { assemble_spades } from "./modules/assembly.nf"
include { assemble_shovill } from "./modules/assembly.nf"
include { assemble_spades_meta } from "./modules/assembly.nf"

include { map_reads } from "./modules/mapping.nf"
include { map_stats } from "./modules/mapping_statistics.nf"
include { sort_mapping } from "./modules/sort_mapping.nf"
include { stats_compile as map_stats_compile } from "./modules/statistics_compilation.nf"
include { stats_compile as assembly_stats_compile } from "./modules/statistics_compilation.nf"
include { stats_compile as classification_stats_compile } from "./modules/statistics_compilation.nf"
include { stats_compile as abundance_S_stats_compile } from "./modules/statistics_compilation.nf"
include { stats_compile as abundance_G_stats_compile } from "./modules/statistics_compilation.nf"

include { map_reads as map_reads_mlst } from "./modules/mapping.nf"

include { variantcall_single } from "./modules/variantcall.nf"
include { variantcall_multi } from "./modules/variantcall.nf"

include { variantcall_single_statistics } from "./modules/variantcall_statistics.nf"
include { compile_variant_counts } from "./modules/variantcall_statistics.nf"

include { rename_skesa } from "./modules/rename_assembly.nf" // actually skesa is not run in assembly, but shovill-spades
include { rename_skesa as rename_spades_meta } from "./modules/rename_assembly.nf"

include { compress_skesa } from "./modules/compress_assembly.nf" // actually skesa is not run in assembly, but shovill-spades
include { compress_skesa as compress_spades_meta } from "./modules/compress_assembly.nf"

include { qc_assembly } from "./modules/qc_assembly.nf"

include { classify } from "./modules/classify.nf"
include { classify as classify_qc } from "./modules/classify.nf"
include { classify_contigs } from "./modules/classify.nf"
include { abundance as abundance_species} from "./modules/classify"
include { abundance as abundance_genus} from "./modules/classify"
include { stats_compile as classification_qc_stats_compile } from "./modules/statistics_compilation.nf"
include { stats_compile as classification_contigs_stats_compile } from "./modules/statistics_compilation.nf"


include { count_reads as count_trimmed_reads } from "./modules/read_count.nf"
include { compile_read_counts as compile_trimmed_read_counts } from "./modules/read_count.nf"
include { count_reads as count_filtered_reads } from "./modules/read_count.nf"
include { compile_read_counts as compile_filtered_read_counts } from "./modules/read_count.nf"

include { filter_reads } from "./modules/filter_reads.nf"
include { filter_contigs } from "./modules/filter_reads.nf"

include { compress_reads } from "./modules/compress_reads.nf"

include { abricate as abricate_argannot} from "./modules/amr_identification.nf"
include { abricate as abricate_plasmidfinder } from "./modules/amr_identification.nf"
include { abricate as abricate_megares } from "./modules/amr_identification.nf"
include { abricate as abricate_ncbi } from "./modules/amr_identification.nf"
include { abricate as abricate_card } from "./modules/amr_identification.nf"
include { abricate as abricate_vfdb } from "./modules/amr_identification.nf"
include { abricate as abricate_resfinder } from "./modules/amr_identification.nf"

include { abricate as abricate_argannot_meta} from "./modules/amr_identification.nf"
include { abricate as abricate_plasmidfinder_meta } from "./modules/amr_identification.nf"
include { abricate as abricate_megares_meta } from "./modules/amr_identification.nf"
include { abricate as abricate_ncbi_meta } from "./modules/amr_identification.nf"
include { abricate as abricate_card_meta } from "./modules/amr_identification.nf"
include { abricate as abricate_vfdb_meta } from "./modules/amr_identification.nf"
include { abricate as abricate_resfinder_meta } from "./modules/amr_identification.nf"



include { amr_extract as amr_extract_argannot} from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_plasmidfinder } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_megares } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_ncbi } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_card } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_vfdb } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_resfinder } from "./modules/amr_extraction.nf"

include { amr_extract as amr_extract_argannot_meta} from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_plasmidfinder_meta } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_megares_meta } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_ncbi_meta } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_card_meta } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_vfdb_meta } from "./modules/amr_extraction.nf"
include { amr_extract as amr_extract_resfinder_meta } from "./modules/amr_extraction.nf"

include { type_contigs_sistr } from "./modules/typing.nf"

include { metagenomics_setup } from "./modules/metagenomics.nf"
include { metagenomics } from "./modules/metagenomics.nf"
 
/************************** 
* WORKFLOW
**************************/

// read in csv file (sep=; due to Excel)
// sample;raw_fwd;raw_rev;raw_folder;trimmed_folder

ch_infiles = Channel
                .fromPath(inputCSV)
                .splitCsv(header: true, sep: ',')
                .map {row -> tuple(
                                    row.sample, 
                                    file(row.raw_folder + "/" + row.raw_fwd),
                                    file(row.raw_folder + "/" + row.raw_rev),
                                    )}


// this workflow takes care of the assembly of isolates and doesn't do any species filtering
workflow ASSEMBY_NOFILTER {

    take:
        ch_infiles

    main:

        // index reference genome
        // TOOL: BWA
        // INPUT: reference genome
        // OUTPUT: trimmed reads
        ch_indexed_reference_genome = index_reference(tuple "reference", reference_fasta)

        // trim reads / run QC
        // INPUT: raw reads
        // OUTPUT: trimmed reads
        // OUTPUT: QC report
        ch_trimmed = trim(ch_infiles) // run trimming

        // Reads mapping
        // TOOL: bwa
        // INPUT: trimmed reads (no species filter)
        // OUTPUT: mapped reads (no species filter)
        ch_mapped_reads = map_reads(ch_indexed_reference_genome, ch_trimmed.fastq, "unfiltered")

        // Read mapping stats
        // TOOL: samtools
        // INPUT: mapped reads (no species filter)
        // OUTPUT: mapping statistics (no species filter)
        ch_mapped_stats = map_stats(ch_mapped_reads.bam, "unfiltered")

        // Read mapping stats compilation
        // TOOL: multiqc
        // INPUT: mapping statistics (no species filter)
        // OUTPUT: mapping statistics report (no species filter)
        ch_mapped_stats_report = map_stats_compile(ch_mapped_stats.collect{ it[1] }, "mapping")
        
        // assemble trimmed reads 
        // TOOL: spades
        // INPUT: trimmed reads
        // OUTPUT: assemblies (no species filter)
        ch_assembled_reads = assemble_shovill(ch_trimmed.fastq, "unfiltered")

        // rename, compress and index assemblies
        // TOOL: -, gzip, bwa
        // INPUT: assemblies (no species filter)
        // OUTPUT: compressed assemblies (no species filter)
        ch_renamed_assembly = rename_skesa(ch_assembled_reads, "unfiltered")
        ch_compressed_assembly = compress_skesa(ch_renamed_assembly.contigs, "unfiltered")
        ch_indexed_assembly = index_assembly(ch_compressed_assembly.contigs)

        // QC assemblies
        // TOOL: Quast
        // INPUT: compressed assemblies (no species filter)
        // OUTPUT: assemblies QC report (no species filter)
        ch_checked_assembly = qc_assembly(ch_compressed_assembly.contigs, reference_fasta, reference_gff, "unfiltered")

        // Variant calling 
        // TOOL: Freebayes
        // INPUT: sorted mapped reads (species filter)
        // OUTPUT 1: single sample variant call (species filter)
        // OUTPUT 2: multi sample variant call (species filter)
        ch_snp_single = variantcall_single( ch_mapped_reads.bam_sorted, 
                                            ch_mapped_reads.bai, 
                                            reference_fasta, 
                                            "unfiltered")
        ch_snp_multi  = variantcall_multi(  ch_mapped_reads.bam_sorted.collect{ it[1] }, 
                                            ch_mapped_reads.bai.collect{ it[1] }, 
                                            reference_fasta, 
                                            "unfiltered")

        // Variant filtering -- not implemented
        // TOOL: bcftools
        // INPUT: variant calls (species filter)
        // OUTPUT: filtered variant calls (species filter)
        
        
        // Classify assemblies
        // TOOL: Kraken
        // INPUT: compressed assemblies (no species filter)
        // OUTPUT: classified assemblies (no species filter)
        if (params.krakendb) {
            ch_classified_assembly = classify_contigs(ch_compressed_assembly.contigs, krakendb, "unfiltered")
        }
        
        // Read assembly stats compilation
        // TOOL: multiqc
        // INPUT: assembly statistics (no species filter)
        // OUTPUT: assembly statistics report (no species filter)
        ch_assembly_stats_report = assembly_stats_compile(ch_checked_assembly.report_tsv.collect{ it[1] }, "assembly")

        // detect AMR genes
        // TOOL: abricate
        // INPUT: compressed assemblies (no species filter), abricate database
        // OUTPUT: listed AMR genes
        ch_amr_genes_argannot = abricate_argannot(ch_renamed_assembly.contigs, "argannot", "unfiltered")
        ch_amr_genes_plasmidfinder = abricate_plasmidfinder(ch_renamed_assembly.contigs, "plasmidfinder", "unfiltered")
        ch_amr_genes_megares = abricate_megares(ch_renamed_assembly.contigs, "megares", "unfiltered")
        ch_amr_genes_ncbi = abricate_ncbi(ch_renamed_assembly.contigs, "ncbi", "unfiltered")
        ch_amr_genes_card = abricate_card(ch_renamed_assembly.contigs, "card", "unfiltered")
        ch_amr_genes_vfdb = abricate_vfdb(ch_renamed_assembly.contigs, "vfdb", "unfiltered")
        ch_amr_genes_resfinder = abricate_resfinder(ch_renamed_assembly.contigs, "resfinder", "unfiltered")


        // extract AMR gene sequences from assemblies
        // TOOL: seqtk
        // INPUT: compressed assemblies + AMR gene positions (start position - 1!)
        // OUTPUT: AMR gene sequences
        ch_amr_genes_sequence_argannot =        amr_extract_argannot(       ch_amr_genes_argannot.join(ch_renamed_assembly.contigs), "argannot", "unfiltered")
        ch_amr_genes_sequence_plasmidfinder =   amr_extract_plasmidfinder(  ch_amr_genes_plasmidfinder.join(ch_renamed_assembly.contigs), "plasmidfinder", "unfiltered")
        ch_amr_genes_sequence_megares =         amr_extract_megares(        ch_amr_genes_megares.join(ch_renamed_assembly.contigs), "megares", "unfiltered")
        ch_amr_genes_sequence_ncbi =            amr_extract_ncbi(           ch_amr_genes_ncbi.join(ch_renamed_assembly.contigs), "ncbi", "unfiltered")
        ch_amr_genes_sequence_card =            amr_extract_card(           ch_amr_genes_card.join(ch_renamed_assembly.contigs), "card", "unfiltered")
        ch_amr_genes_sequence_vfdb =            amr_extract_vfdb(           ch_amr_genes_vfdb.join(ch_renamed_assembly.contigs), "vfdb", "unfiltered")
        ch_amr_genes_sequence_resfinder =       amr_extract_resfinder(      ch_amr_genes_resfinder.join(ch_renamed_assembly.contigs), "resfinder", "unfiltered")

        // run sample typing
        // TOOL: sistr
        // INPUT: contigs (uncompressed, no species filter)
        // OUPUT: typing reports (no species filter)
        if (params.typing) {
            ch_typed_assembly = type_contigs_sistr(ch_renamed_assembly.contigs.filter{ it[1].size() > 0 }, "unfiltered")
        }

}

workflow ASSEMBLY_FILTER {

    take:
        ch_infiles

    main:

        // index reference genome
        // TOOL: BWA
        // INPUT: reference genome
        // OUTPUT: indexed genome
        ch_indexed_reference_genome = index_reference(tuple "reference", reference_fasta)

        // trim reads / run QC
        // INPUT: raw reads
        // OUTPUT: trimmed reads
        // OUTPUT: QC report
        ch_trimmed = trim(ch_infiles) // run trimming

        // count trimmed reads 
        // TOOL: awk
        // INPUT: trimmed reads
        // OUTPUT: trimmed read counts
        ch_counted_trimmed_reads = count_trimmed_reads(ch_trimmed.fastq, "01_trimming")
        ch_compiled_counted_trimmed_reads = compile_trimmed_read_counts(ch_counted_trimmed_reads.read_count.collect{ it[1] }, "trimmed")

        // classify reads 
        // TOOL: Kraken2
        // INPUT: trimmed reads
        // OUTPUT: classified reads
        ch_classified = classify(ch_trimmed.fastq, krakendb, "filtered")

        // filter classify reads 
        // TOOL: krakentools
        // INPUT: classified reads
        // OUTPUT: filtered & classified reads
        ch_filtered_reads = filter_reads(ch_classified.fastq, ch_classified.kraken_output, ch_classified.kraken_report, taxid)

        // QC filter classify reads 
        // TOOL: kraken
        // INPUT: classified + filtered reads
        // OUTPUT: QC report of filtered reads
        ch_filtered_reads_qc = classify_qc(ch_filtered_reads.fastq, krakendb_contigs, "qc")

        // classify reads statistics
        // TOOL: multiqc
        // INPUT: classified reads & reports
        // OUTPUT: report on classification results
        ch_filtered_reads_statistics = classification_stats_compile(ch_classified.kraken_report.collect{ it[1] }, "classification")

        // classify reads QC statistics
        // TOOL: multiqc
        // INPUT: QC of classified reads & reports
        // OUTPUT: report on QC classification results
        ch_filtered_reads_statistics = classification_qc_stats_compile(ch_filtered_reads_qc.kraken_report.collect{ it[1] }, "classification_QC")


        // abundance estimation
        // TOOL: Bracken
        // INPUT: classification report
        // OUTPUT: abundance estimation
        ch_abundance_species = abundance_species(ch_classified.kraken_report, krakendb, "S")
        ch_abundance_genus   = abundance_genus(ch_classified.kraken_report, krakendb, "G")

        // abundance statistics
        // TOOL: multiqc
        // INPUT: abundance estimation
        // OUTPUT: report on abundance estimation
        ch_abundance_species_report = abundance_S_stats_compile(ch_abundance_species.bracken_report.collect{ it[1] }, "abundance_species")
        ch_abundance_genus_report   = abundance_G_stats_compile(ch_abundance_genus.bracken_report.collect{ it[1] }, "abundance_genus")

        // compress filtered & classify reads 
        // TOOL: gzip
        // INPUT: filtered & classified reads
        // OUTPUT: compressed & filtered & classified reads
        ch_compressed_reads = compress_reads(ch_filtered_reads.fastq)

        // count filtered & classify reads 
        // TOOL: awk
        // INPUT: filtered & classified reads
        // OUTPUT: filtered read counts
        ch_counted_filtered_reads = count_filtered_reads(ch_compressed_reads.fastq, "02_read_filter")
        ch_compiled_counted_filtered_reads = compile_filtered_read_counts(ch_counted_filtered_reads.read_count.collect{ it[1] }, "filtered")


        // Reads mapping
        // TOOL: bwa
        // INPUT: compressed reads (species filter)
        // OUTPUT: mapped reads (species filter)
        ch_mapped_reads = map_reads(ch_indexed_reference_genome, ch_compressed_reads.fastq.filter{ it[1].size() > 0 }, "filtered")

        // Read mapping stats
        // TOOL: samtools
        // INPUT: mapped reads (species filter)
        // OUTPUT: mapping statistics (species filter)
        ch_mapped_stats = map_stats(ch_mapped_reads.bam, "filtered")

        // Read mapping stats compilation
        // TOOL: multiqc
        // INPUT: mapping statistics (species filter)
        // OUTPUT: mapping statistics report (species filter)
        ch_mapped_stats_report = map_stats_compile(ch_mapped_stats.collect{ it[1] }, "mapping")

        // Variant calling 
        // TOOL: Freebayes
        // INPUT: sorted mapped reads (species filter)
        // OUTPUT 1: single sample variant call (species filter)
        // OUTPUT 2: multi sample variant call (species filter)
        ch_snp_single = variantcall_single( ch_mapped_reads.bam_sorted, 
                                            ch_mapped_reads.bai, 
                                            reference_fasta, 
                                            "filtered")
        ch_snp_multi  = variantcall_multi(  ch_mapped_reads.bam_sorted.collect{ it[1] }, 
                                            ch_mapped_reads.bai.collect{ it[1] }, 
                                            reference_fasta, 
                                            "filtered")

        // Variant stats
        // TOOL: bcftools
        // INPUT: variant calls (species filter)
        // OUTPUT: variant statistics (species filter)
        ch_variant_stats = variantcall_single_statistics(ch_snp_single.variants)
        ch_compiled_variant_stats = compile_variant_counts(ch_variant_stats.collect{ it[1] })

        // Variant filtering -- not implemented
        // TOOL: bcftools
        // INPUT: variant calls (species filter)
        // OUTPUT: filtered variant calls (species filter)

        if (params.mlst_mapping) {
            // index mlst reference 
            // TOOL: BWA
            // INPUT: mlst reference
            // OUTPUT: indexed mlst reference
            ch_indexed_mlst = index_mlst(tuple "mlst", params.mlst_mapping)

            // Reads mapping
            // TOOL: bwa
            // INPUT: compressed reads (species filter)
            // OUTPUT: mapped reads (species filter)
            ch_mapped_reads2mlst = map_reads_mlst(ch_indexed_mlst, ch_compressed_reads.fastq.filter{ it[1].size() > 0 }, "filtered_mlst")
        }
        
        // assemble trimmed reads 
        // TOOL: spades
        // INPUT: trimmed filtered reads
        // OUTPUT: assemblies (species filter)
        ch_assembled_reads = assemble_shovill(ch_compressed_reads.fastq.filter{ it[1].size() > 0 }, "filtered")

        // rename, compress and index assemblies
        // TOOL: -, gzip, bwa
        // INPUT: assemblies (species filter)
        // OUTPUT: compressed assemblies (species filter)
        ch_renamed_assembly = rename_skesa(ch_assembled_reads, "filtered")
        ch_compressed_assembly = compress_skesa(ch_renamed_assembly.contigs, "filtered")
        ch_indexed_assembly = index_assembly(ch_compressed_assembly.contigs)

        // QC assemblies
        // TOOL: Quast
        // INPUT: compressed assemblies (species filter)
        // OUTPUT: assemblies QC report (species filter)
        ch_checked_assembly = qc_assembly(ch_compressed_assembly.contigs, reference_fasta, reference_gff, "filtered")


        // QC assemblies / taxonomy check
        // TOOL: Kraken
        // INPUT: compressed assemblies (species filter)
        // OUTPUT: taxonomic classification of contigs (species filter) 
        ch_classified_assembly = classify_contigs(ch_renamed_assembly.contigs.filter{ it[1].size() > 0 }, krakendb_contigs, "filtered")

        // classify reads QC statistics
        // TOOL: multiqc
        // INPUT: QC of classified reads & reports
        // OUTPUT: report on QC classification results
        ch_filtered_reads_statistics = classification_contigs_stats_compile(ch_classified_assembly.kraken_report.collect{ it[1] }, "classification_contigs")

        // Read assembly stats compilation
        // TOOL: multiqc
        // INPUT: assembly statistics (species filter)
        // OUTPUT: assembly statistics report (species filter)
        ch_assembly_stats_report = assembly_stats_compile(ch_checked_assembly.report_tsv.collect{ it[1] }, "assembly")

        // filter classify contigs 
        // TOOL: krakentools
        // INPUT: classified contigs
        // OUTPUT: filtered & classified contigs
        ch_filtered_reads = filter_contigs(ch_classified_assembly.fasta, ch_classified_assembly.kraken_output, ch_classified_assembly.kraken_report, taxid, "filtered")

        // detect AMR genes
        // TOOL: abricate
        // INPUT: compressed assemblies (no species filter), abricate database
        // OUTPUT: listed AMR genes
        ch_amr_genes_argannot       = abricate_argannot(ch_renamed_assembly.contigs, "argannot", "filtered")
        ch_amr_genes_plasmidfinder  = abricate_plasmidfinder(ch_renamed_assembly.contigs, "plasmidfinder", "filtered")
        ch_amr_genes_megares        = abricate_megares(ch_renamed_assembly.contigs, "megares", "filtered")
        ch_amr_genes_ncbi           = abricate_ncbi(ch_renamed_assembly.contigs, "ncbi", "filtered")
        ch_amr_genes_card           = abricate_card(ch_renamed_assembly.contigs, "card", "filtered")
        ch_amr_genes_vfdb           = abricate_vfdb(ch_renamed_assembly.contigs, "vfdb", "filtered")
        ch_amr_genes_resfinder      = abricate_resfinder(ch_renamed_assembly.contigs, "resfinder", "filtered")

        // extract AMR gene sequences from assemblies
        // TOOL: seqtk
        // INPUT: compressed assemblies + AMR gene positions (start position - 1!)
        // OUTPUT: AMR gene sequences
        ch_amr_genes_sequence_argannot =        amr_extract_argannot(       ch_amr_genes_argannot.join(ch_renamed_assembly.contigs), "argannot", "filtered")
        ch_amr_genes_sequence_plasmidfinder =   amr_extract_plasmidfinder(  ch_amr_genes_plasmidfinder.join(ch_renamed_assembly.contigs), "plasmidfinder", "filtered")
        ch_amr_genes_sequence_megares =         amr_extract_megares(        ch_amr_genes_megares.join(ch_renamed_assembly.contigs), "megares", "filtered")
        ch_amr_genes_sequence_ncbi =            amr_extract_ncbi(           ch_amr_genes_ncbi.join(ch_renamed_assembly.contigs), "ncbi", "filtered")
        ch_amr_genes_sequence_card =            amr_extract_card(           ch_amr_genes_card.join(ch_renamed_assembly.contigs), "card", "filtered")
        ch_amr_genes_sequence_vfdb =            amr_extract_vfdb(           ch_amr_genes_vfdb.join(ch_renamed_assembly.contigs), "vfdb", "filtered")
        ch_amr_genes_sequence_resfinder =       amr_extract_resfinder(      ch_amr_genes_resfinder.join(ch_renamed_assembly.contigs), "resfinder", "filtered")



        // run sample typing
        // TOOL: sistr
        // INPUT: contigs (uncompressed, no species filter)
        // OUPUT: typing reports (no species filter)
        if (params.typing) {
            ch_typed_assembly = type_contigs_sistr(ch_renamed_assembly.contigs.filter{ it[1].size() > 0 }, "filtered")
        }

        // assemble trimmed reads 
        // TOOL: spades -meta
        // INPUT: trimmed filtered reads
        // OUTPUT: assemblies meta
        ch_assembled_meta_reads = assemble_spades_meta(ch_trimmed.fastq.filter{ it[1].size() > 0 }, "meta", "--meta")

        // rename, compress and index assemblies
        // TOOL: -, gzip, bwa
        // INPUT: assemblies (species filter)
        // OUTPUT: compressed assemblies (species filter)
        ch_renamed_assembly_meta    = rename_spades_meta(ch_assembled_meta_reads.contigs, ch_assembled_meta_reads.graph, "meta")
        ch_compressed_assembly_meta = compress_spades_meta(ch_renamed_assembly_meta.contigs, "meta")
        ch_indexed_assembly_meta    = index_assembly_meta(ch_compressed_assembly_meta.contigs)

        // detect AMR genes in metagenomic assemblies
        // TOOL: abricate
        // INPUT: compressed assemblies (metagenomes), abricate database
        // OUTPUT: listed AMR genes
        ch_amr_genes_argannot_meta =         abricate_argannot_meta(ch_renamed_assembly_meta.contigs, "argannot", "meta")
        ch_amr_genes_plasmidfinder_meta =    abricate_plasmidfinder_meta(ch_renamed_assembly_meta.contigs, "plasmidfinder", "meta")
        ch_amr_genes_megares_meta =          abricate_megares_meta(ch_renamed_assembly_meta.contigs, "megares", "meta")
        ch_amr_genes_ncbi_meta =             abricate_ncbi_meta(ch_renamed_assembly_meta.contigs, "ncbi", "meta")
        ch_amr_genes_card_meta =             abricate_card_meta(ch_renamed_assembly_meta.contigs, "card", "meta")
        ch_amr_genes_vfdb_meta =             abricate_vfdb_meta(ch_renamed_assembly_meta.contigs, "vfdb", "meta")
        ch_amr_genes_resfinder_meta =        abricate_resfinder_meta(ch_renamed_assembly_meta.contigs, "resfinder", "meta")

        // extract AMR gene sequences from metagenome assemblies
        // TOOL: seqtk
        // INPUT: compressed metagenome assemblies + AMR gene positions (start position - 1!)
        // OUTPUT: AMR gene sequences (metagenome)
        ch_amr_genes_sequence_argannot_meta =        amr_extract_argannot_meta(       ch_amr_genes_argannot_meta.join(ch_renamed_assembly_meta.contigs), "argannot", "meta")
        ch_amr_genes_sequence_plasmidfinder_meta =   amr_extract_plasmidfinder_meta(  ch_amr_genes_plasmidfinder_meta.join(ch_renamed_assembly_meta.contigs), "plasmidfinder", "meta")
        ch_amr_genes_sequence_megares_meta =         amr_extract_megares_meta(        ch_amr_genes_megares_meta.join(ch_renamed_assembly_meta.contigs), "megares", "meta")
        ch_amr_genes_sequence_ncbi_meta =            amr_extract_ncbi_meta(           ch_amr_genes_ncbi_meta.join(ch_renamed_assembly_meta.contigs), "ncbi", "meta")
        ch_amr_genes_sequence_card_meta =            amr_extract_card_meta(           ch_amr_genes_card_meta.join(ch_renamed_assembly_meta.contigs), "card", "meta")
        ch_amr_genes_sequence_vfdb_meta =            amr_extract_vfdb_meta(           ch_amr_genes_vfdb_meta.join(ch_renamed_assembly_meta.contigs), "vfdb", "meta")
        ch_amr_genes_sequence_resfinder_meta =       amr_extract_resfinder_meta(      ch_amr_genes_resfinder_meta.join(ch_renamed_assembly_meta.contigs), "resfinder", "meta")



        if (params.metaphlan_db) {
            // ch_metagenome_setup = metagenomics_setup() # needs whitelist entry
            ch_metagenome = metagenomics(ch_trimmed.fastq, "filtered", params.metaphlan_db)
        }
        // // All stats compilation
        // // TOOL: multiqc
        // // INPUT: trimming statistics (species filter)
        // // INPUT: classification statistics (species filter)
        // // INPUT: mapping statistics (species filter)
        // // INPUT: assembly statistics (species filter)
        // // OUTPUT: all statistics report (species filter)
        // ch_mapped_stats_report = all_stats_compile( ch_trimmed.json.collect{ it[1] },
        //                                             ch_classified.kraken_report.collect{ it[1] },
        //                                             ch_mapped_stats.collect{ it[1] }, 
        //                                             ch_checked_assembly.report_tsv.collect{ it[1] },
        //                                             "all")

}

workflow {

    if (params.krakendb && params.taxid) {
        print("chose to assemble kraken filtered reads")
        ASSEMBLY_FILTER(ch_infiles)
    }
    else {
        print("chose to assemble unfiltered reads")
        ASSEMBY_NOFILTER(ch_infiles)
    }
    

}

// this process may also be moved into a module, but is located here to show how "inline" processes are defined.
process trim {

    label "fastp"

    conda "${projectDir}/envs/fastp.yaml"

    // this line will publish all "fastq.gz" files to the output directory
    // publishDir "${params.analysesdir}/01_trimming/${sample}", pattern: "*.fastq.gz", mode: "copy"

    // if you want to publish all output that was defined in the "output:" section further down, one leaves out the "pattern:"
    // we're using this one, because we also want to publish the fastp trimming reports
    publishDir "${params.analysesdir}/01_trimming/${sample}", mode: "copy"

    input:
        tuple val(sample), path(read1), path(read2)

    // each output of fastp is tagged with an "emit" to directly reference it for downstream usage
    output:
        tuple val(sample), path("${sample}.R*.trimmed.fastq.gz"), emit: fastq
        tuple val(sample), path("${sample}.trim_report.json"), emit: json
        tuple val(sample), path("${sample}.trim_report.html"), emit: html

    // SLURM cluster options
    cpus 10
    memory "3 GB"
    time "1h"
    tag "trim_${sample}"

    // here we write what is going to be executed
    script:
        """
            fastp \
                --in1 ${read1} \
                --in2 ${read2} \
                --out1 ${sample}.R1.trimmed.fastq.gz \
                --out2 ${sample}.R2.trimmed.fastq.gz \
                --trim_poly_g \
                --json ${sample}.trim_report.json \
                --html ${sample}.trim_report.html \
                --thread ${task.cpus} \
                --report_title ${sample}
        """

    // we can start the pipeline in "--stub" mode. Then this stub is going to be executed that basically checks if the workflow processes fit together
    stub:
        """
            touch ${sample}.R{1,2}.trimmed.fastq.gz ${sample}.trim_report.json ${sample}.trim_report.html
        """

}
