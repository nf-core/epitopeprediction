//
// Check input samplesheet and get read channels
//

include { SYFPEITHI } from '../../modules/local/syfpeithi'
include { MHCFLURRY } from '../../modules/local/mhcflurry'
include { MHCNUGGETS } from '../../modules/local/mhcnuggets'
include { NETMHCPAN } from '../../modules/local/netmhcpan'
include { NETMHCIIPAN } from '../../modules/local/netmhciipan'
include { EXTERNAL_TOOLS_IMPORT} from '../../modules/local/external_tools_import'
include { MERGE_PREDICTIONS } from '../../modules/local/merge_predictions'
include { PREPARE_PREDICTION_INPUT } from '../../modules/local/prepare_prediction_input'

workflow MHC_BINDING_PREDICTION {
    take:
        ch_peptides

    main:
        ch_versions = Channel.empty()
        ch_combined_predictions = Channel.empty()

        if (params.tools.isEmpty()) { exit 1, "No valid tools specified." }

        tools = params.tools.split(',')

        // Create a channel element for each tool-tsv combination and save tool as key in meta map
        ch_peptides
            .flatMap { meta, file ->
                tools.collect { tool ->
                    def newMeta = meta.clone()
                    newMeta.put('tool', tool)
                    tuple(newMeta, file)
                }
            }
            .set{ ch_peptides_to_prepare }

        //prepare the input file
        PREPARE_PREDICTION_INPUT( ch_peptides_to_prepare )
            .prepared
            .branch {
                meta, peptide_file ->
                    syfpeithi : meta.tool == 'syfpeithi'
                        return [meta, peptide_file]
                    mhcflurry : meta.tool == 'mhcflurry'
                        return [meta, peptide_file]
                    mhcnuggets : meta.tool == 'mhcnuggets'
                        return [meta, peptide_file]
                    netmhcpan: meta.tool == 'netmhcpan'
                        return [meta, peptide_file]
                    netmhciipan: meta.tool == 'netmhciipan'
                        return [meta, peptide_file]
                    }
            .set{ ch_prediction_input }

        SYFPEITHI ( ch_prediction_input.syfpeithi )
        ch_versions = ch_versions.mix(SYFPEITHI.out.versions)

        MHCFLURRY ( ch_prediction_input.mhcflurry )
        ch_versions = ch_versions.mix(MHCFLURRY.out.versions)

        MHCNUGGETS ( ch_prediction_input.mhcnuggets )
        ch_versions = ch_versions.mix(MHCNUGGETS.out.versions)

        if ( "netmhc" in tools )
        {
            ch_netmhc_tool = EXTERNAL_TOOLS_IMPORT( parse_netmhc_params("netmhcpan", "4.1") )
        }
        else
        {
            ch_netmhc_tool = Channel.empty()
        }
        NETMHCPAN ( ch_prediction_input.netmhcpan.combine(ch_netmhc_tool) )
        ch_versions = ch_versions.mix(NETMHCPAN.out.versions)

        // TODO: External tools import for netmhciipan
        NETMHCIIPAN ( ch_prediction_input.netmhciipan.combine(ch_netmhc_tool) )
        ch_versions = ch_versions.mix(NETMHCIIPAN.out.versions)

    ch_combined_predictions
        .join( SYFPEITHI.out.predicted, remainder: true )
        .join( MHCFLURRY.out.predicted, remainder: true )
        .join( MHCNUGGETS.out.predicted, remainder: true )
        .join( NETMHCPAN.out.predicted, remainder: true )
        .join( NETMHCIIPAN.out.predicted, remainder: true )
        .map { it.findAll { item -> item != null } }
        .map { meta, predicted_file -> [meta.subMap(['sample', 'alleles', 'mhcclass', 'inputtype']), predicted_file] }
        .groupTuple()
        .join( ch_peptides )
        .set{ ch_combined_predictions_meta }

    //merge the prediction output of all tools into one output merged_prediction.tsv
    MERGE_PREDICTIONS (ch_combined_predictions_meta)
    ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    emit:
    //predicted = MERGE_PREDICTIONS.out.merged
    predicted = Channel.empty()
    versions = ch_versions
}


// Functions
def parse_netmhc_params(tool_name, tool_version) {
    // Check if the _path parameter was set for this tool
    if (!params["${tool_name}_path"])
    {
        error("--${tool_name}_path not specified, but --tools contains ${tool_name}. Both have to be specified to enable ${tool_name}. Ignoring.")
    }
    else if (params["${tool_name}_path"])
    {
    // Import mandatory netmhc metadata
    def jsonSlurper = new groovy.json.JsonSlurper()
    def external_tools_meta = jsonSlurper.parse(file(params.external_tools_meta, checkIfExists: true))
    def entry = external_tools_meta[tool_name][tool_version]

    if (params["netmhc_system"] == 'darwin') {
        entry = external_tools_meta["${tool_name}_darwin"][tool_version]
    }
    // If so, add the tool name and user installation path to the external tools import channel
    ch_nonfree_paths = Channel.empty()
    ch_nonfree_paths.bind([
        tool_name,
        entry.version,
        entry.software_md5,
        file(params["${tool_name}_path"], checkIfExists:true),
        file(entry.data_url),
        entry.data_md5,
        entry.binary_name
    ])

    return ch_nonfree_paths
    }
}
