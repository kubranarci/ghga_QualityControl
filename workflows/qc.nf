/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { TABIX_TABIX            } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS         } from '../modules/nf-core/bcftools/stats/main'
include { RSEQC_BAMSTAT          } from '../modules/nf-core/rseqc/bamstat/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { NGSBITS_SAMPLEGENDER   } from '../modules/nf-core/ngsbits/samplegender/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    // create reference channels ////

    ch_fasta       = Channel.fromPath(params.fasta, checkIfExists: true)
                        .map{ fasta -> tuple([id: fasta.getSimpleName()], fasta) }.collect()


    // create empty channels for versions and multiqc files
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet.branch{ it ->
        def meta = it[0]
        step1: meta.step == 1
        step2: meta.step == 2
        step3: meta.step == 3
        other:false
    }.set{samplesheet}

    samplesheet.step1.view()
    // Runs FASTQC if samplesheet has fastq files
    FASTQC (
        samplesheet.step1
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Runs FASTP if samplesheet has fastq files
    save_trimmed_fail = false
    save_merged = false
    FASTP (
        samplesheet.step1,
        [],
        false,
        save_trimmed_fail,
        save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions)

    if (params.method in ["wgs", "wes", "tes"]){
        // Runs MOSDEPTH if samplesheet has bam/cram files (targetted, wes, wgs will be different)
        MOSDEPTH(
            samplesheet.step2.map{meta, file, index -> tuple(meta, file, index, []) },
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

        // Runs SAMTOOLS_STATS if samplesheet has bam/cram files
        SAMTOOLS_STATS(
            samplesheet.step2,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    }

    if (params.method.contains("rna")){
        RSEQC_BAMSTAT(
            samplesheet.step2
        )
        ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)

    }

    if (params.predict_sex && params.method.contains("wgs") ){

        // this can be replaced if more than once tool will use fai
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[],[]],
            false
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)


        // Predict sex of the samples if not given
        NGSBITS_SAMPLEGENDER(
            samplesheet.step2,
            ch_fasta,
            SAMTOOLS_FAIDX.out.fai,
            params.samplegender_method ?: 'xy'
        )
        ch_multiqc_files = ch_multiqc_files.mix(NGSBITS_SAMPLEGENDER.out.tsv.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions)
    }


    // Runs BCFTOOLS_STATS if samplesheet has vcf files (targetted, wes, wgs will be different)
    // TODO add support for vcf and bcf-bcf.gz
    TABIX_TABIX(
        samplesheet.step3
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    regions = [[],[]]
    targets = [[],[]]
    samples = [[],[]]
    exons   = [[],[]]

    BCFTOOLS_STATS(
        samplesheet.step3.join(TABIX_TABIX.out.tbi),
        regions,
        targets,
        samples,
        exons,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'qc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
