/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_magdown_pipeline'
include { ABRICATE_RUN           } from '../modules/nf-core/abricate/run/main'
include { ABRICATE_SUMMARY       } from '../modules/nf-core/abricate/summary/main'
include { PLASMIDFINDER          } from '../modules/nf-core/plasmidfinder/main'
include { PLASMIDFINDER_SUMMARY  } from '../modules/local/plasmidfinder_summary/main'
include { PLASMIDFINDER_MATRIX   } from '../modules/local/plasmidfinder_matrix/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAGDOWN {

    take:
    ch_samplesheet // channel: [ val(meta), path(magdir) ] from samplesheet
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Expand each sample's MAGdir into one channel entry per FASTA file
    // Emits: [ val(meta + [mag_id: <filename>]), path(fasta) ]
    //
    def ch_fastas = ch_samplesheet
        .flatMap { meta, magdir ->
            magdir.listFiles()
                .findAll { it.name =~ /\.f(a|asta|na)(\.gz)?$/ }
                .collect { fasta ->
                    def mag_id = fasta.name.replaceAll(/\.f(a|asta|na)(\.gz)?$/, '')
                    [ meta + [mag_id: mag_id, id: mag_id, sample: meta.id], fasta ]
                }
        }

    //
    // MODULE: Run Abricate on each MAG
    //
    ABRICATE_RUN(ch_fastas, [])

    //ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions_abricate.map { _process, _tool, version -> version })

    //
    // Group per-MAG reports back by sample for summary
    //
    def ch_abricate_reports_by_sample = ABRICATE_RUN.out.report
        .map { meta, report -> [ [id: meta.sample], report ] }
        .groupTuple()

    ABRICATE_SUMMARY(ch_abricate_reports_by_sample)

    // 
    // MODULE: Run plasmidfinder on each MAG
    //

    PLASMIDFINDER(ch_fastas)

    //
    // Group per-MAG plasmidfinder outputs back by sample for summary
    //
    def ch_plasmidfinder_tsvs_by_sample = PLASMIDFINDER.out.tsv
        .map { meta, tsv -> [ [id: meta.sample], tsv ] }
        .groupTuple()

    PLASMIDFINDER_SUMMARY(ch_plasmidfinder_tsvs_by_sample)

    //
    // Collect per-sample summaries into one samples x plasmids matrix
    //
    def ch_plasmidfinder_summaries = PLASMIDFINDER_SUMMARY.out.summary
        .map { _meta, summary -> summary }
        .collect()
        .map { summaries -> [ [id: 'all_samples'], summaries ] }

    PLASMIDFINDER_MATRIX(ch_plasmidfinder_summaries)

    //ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions_plasmidfinder.map { _process, _tool, version -> version })


    // 
    // MODULE: MULTIQC
    //

    ch_multiqc_files = ch_multiqc_files.mix(ABRICATE_RUN.out.report.map { _meta, file -> file })

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'magdown_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'magdown'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                                    // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
