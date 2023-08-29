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


workflow MHC_BINDING_PREDICTION {
    take:
        metadata_and_file

    main:
        ch_versions = Channel.empty()
        ch_combined_predictions = Channel.empty()

        if (params.tools.isEmpty()) { exit 1, "No valid tools specified." }

        tools = params.tools.split(',')

        if ( "syfpeithi" in tools )
        {
        SYFPEITHI ( metadata_and_file )
        ch_versions = ch_versions.mix(SYFPEITHI.out.versions)
        ch_combined_predictions = ch_combined_predictions.join(SYFPEITHI.out.predicted, remainder: true)
        }
        if ( "mhcflurry" in tools )
        {
        MHCFLURRY ( metadata_and_file )
        ch_versions = ch_versions.mix(MHCFLURRY.out.versions)
        ch_combined_predictions = ch_combined_predictions.join(MHCFLURRY.out.predicted, remainder: true)
        }
        if ( "mhcnuggets" in tools )
        {
        MHCNUGGETS ( metadata_and_file )
        ch_versions = ch_versions.mix(MHCNUGGETS.out.versions)
        ch_combined_predictions = ch_combined_predictions.join(MHCNUGGETS.out.predicted, remainder: true)
        }
        if ( "netmhcpan" in tools )
        {
        EXTERNAL_TOOLS_IMPORT (parse_netmhc_params("netmhcpan", "4.1"))
        NETMHCPAN (metadata_and_file.combine(EXTERNAL_TOOLS_IMPORT.out.nonfree_tools))
        ch_versions = ch_versions.mix(NETMHCPAN.out.versions)
        ch_combined_predictions = ch_combined_predictions.join(NETMHCPAN.out.predicted, remainder: true)
        }
        if ( "netmhciipan" in tools )
        {
        // TODO: External tools import for netmhciipan
        NETMHCIIPAN (metadata_and_file)
        ch_versions = ch_versions.mix(NETMHCIIPAN.out.versions)
        ch_combined_predictions = ch_combined_predictions.join(NETMHCIIPAN.out.predicted, remainder: true)
        }

    //remove the null (it[1]) in the channel output and join metadata and input channel with metadata and output channel
    ch_combined_predictions = ch_combined_predictions.map{ it -> [it[0], it[2..-1]]}.join(metadata_and_file, remainder: true)
    ch_combined_predictions.view()

    //merge the prediction output of all tools into one output merged_prediction.tsv
    MERGE_PREDICTIONS (ch_combined_predictions)

    ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    emit:
    predicted = MERGE_PREDICTIONS.out.merged
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
