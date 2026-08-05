//==============================================================================
// Predict MHC binding for a set of peptides using different predictors.
//==============================================================================

include { PREPARE_PREDICTION_INPUT                   } from '../../../modules/local/prepare_prediction_input'
include { MHCFLURRY                                  } from '../../../modules/local/mhcflurry'
include { MHCNUGGETS;
        MHCNUGGETS as MHCNUGGETSII                   } from '../../../modules/local/mhcnuggets'
include { NETMHCPAN                                  } from '../../../modules/local/netmhcpan'
include { NETMHCIIPAN                                } from '../../../modules/local/netmhciipan'
include { UNPACK_NETMHC_SOFTWARE as NETMHCPAN_IMPORT;
        UNPACK_NETMHC_SOFTWARE as NETMHCIIPAN_IMPORT } from '../../../modules/local/unpack_netmhc_software'
include { MERGE_PREDICTIONS                          } from '../../../modules/local/merge_predictions'

// Input:
//     ch_peptides: Channel of peptides to predict
//     tools: Comma-separated list of tools to use for prediction
//     supported_alleles_json: JSON file with supported alleles for each predictor
//     netmhc_software_meta_json: JSON file with metadata for NetMHC software
// Output:
//     predicted: Channel of predicted MHC binding
//     versions: Channel of software versions

workflow MHC_BINDING_PREDICTION {
    take:
        ch_peptides
        tools
        supported_alleles_json
        netmhc_software_meta_json

    main:
        ch_versions = channel.empty()
        ch_binding_predictors_out = channel.empty()

        validate_tools_param(tools)

        // SPLIT_PEPTIDES already names its chunks <sample>[_<split>]_c<N>, which is unique across the run
        ch_peptides
            .map { meta, file -> [meta + [file_id: file.baseName], file] }
            .set { ch_peptides_to_predict }

        // Fan out one tuple per (tool, chunk) entry from the emitted JSON manifest.
        PREPARE_PREDICTION_INPUT( ch_peptides_to_predict, supported_alleles_json)
            .prepared
            .flatMap { meta, tool_chunks, files ->
                def files_by_name = files.collectEntries { f -> [(f.name): f] }
                parseJson(tool_chunks).collect { entry ->
                    [meta + [tool: entry.tool,
                             alleles_supported: entry.alleles,
                             source_file_id: meta.file_id,
                             file_id: entry.chunk_id ? "${meta.file_id}_${entry.chunk_id}" : meta.file_id],
                     entry.alleles_input,
                     files_by_name[entry.filename]]
                }
            }
            .branch { meta, _alleles_input, _file ->
                mhcflurry    : meta.tool == 'mhcflurry'
                mhcnuggets   : meta.tool == 'mhcnuggets'
                mhcnuggetsii : meta.tool == 'mhcnuggetsii'
                netmhcpan    : meta.tool == 'netmhcpan'
                netmhciipan  : meta.tool == 'netmhciipan'
            }
            .set{ ch_prediction_input }

        // MHCflurry encodes alleles inline in its CSV input, so it doesn't need alleles_input.
        MHCFLURRY ( ch_prediction_input.mhcflurry.map { meta, _alleles_input, file -> [meta, file] } )
        ch_versions = ch_versions.mix(MHCFLURRY.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCFLURRY.out.predicted)

        MHCNUGGETS ( ch_prediction_input.mhcnuggets )
        ch_versions = ch_versions.mix(MHCNUGGETS.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCNUGGETS.out.predicted)

        MHCNUGGETSII ( ch_prediction_input.mhcnuggetsii )
        ch_versions = ch_versions.mix(MHCNUGGETSII.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCNUGGETSII.out.predicted)

        if ( "netmhcpan" in tools.tokenize(",") )
        {
            NETMHCPAN_IMPORT( parse_netmhc_params("netmhcpan", netmhc_software_meta_json) )
            NETMHCPAN ( ch_prediction_input.netmhcpan.combine(NETMHCPAN_IMPORT.out.nonfree_tools) )
            ch_versions = ch_versions.mix(NETMHCPAN.out.versions)
            ch_binding_predictors_out = ch_binding_predictors_out.mix(NETMHCPAN.out.predicted)
        }

        if ( "netmhciipan" in tools.tokenize(",") )
        {
            NETMHCIIPAN_IMPORT( parse_netmhc_params("netmhciipan", netmhc_software_meta_json) )
            NETMHCIIPAN ( ch_prediction_input.netmhciipan.combine(NETMHCIIPAN_IMPORT.out.nonfree_tools) )
            ch_versions = ch_versions.mix(NETMHCIIPAN.out.versions)
            ch_binding_predictors_out = ch_binding_predictors_out.mix(NETMHCIIPAN.out.predicted)
        }

    // Regroup predictor outputs by the original (pre-chunk) peptide file, one MERGE task per source.
    ch_binding_predictors_out
        .map { meta, file ->
            def regroup_meta = meta.findAll { k, _v -> !(k in ['alleles_supported', 'tool', 'source_file_id']) } + [
                file_id: meta.source_file_id ?: meta.file_id,
            ]
            [regroup_meta, file, meta.alleles_supported]
        }
        .groupTuple()                   // → [meta, [files], [alleles_per_file]]
        .join( ch_peptides_to_predict ) // → [meta, [files], [alleles_per_file], source_file]
        .set { ch_binding_predictors_out_meta}

    // Merge predictions from different predictors
    MERGE_PREDICTIONS( ch_binding_predictors_out_meta )
    ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    emit:
    predicted = MERGE_PREDICTIONS.out.merged
    versions = ch_versions
}

//==============================================================================
//                            Auxiliar Functions
//==============================================================================

// Check if supported tools are specified
def validate_tools_param(tools) {
    def valid_tools = [ 'mhcnuggets', 'mhcnuggetsii', 'mhcflurry', 'netmhcpan', 'netmhciipan' ]
    def tool_list = tools.tokenize(',')
    // Validate each tool in tools if it's in valid_tools
    def invalid_tools = tool_list.findAll { tool -> tool.trim() !in valid_tools }
    if (invalid_tools) {
        throw new IllegalArgumentException("Invalid tools found: ${invalid_tools.join(',')}.\nValid tools: ${valid_tools.join(',')}")
    }
}

// Prepare import of NetMHC software
def parse_netmhc_params(tool_name, netmhc_software_meta) {
    // Check if the _path parameter was set for this tool
    if (!params["${tool_name}_path"])
    {
        error("--${tool_name}_path not specified, but --tools contains ${tool_name}. Both have to be specified to enable ${tool_name}. Ignoring.")
    }
    // Import mandatory netmhc metadata
    def jsonSlurper = new groovy.json.JsonSlurper()
    def netmhc_software_meta_map = jsonSlurper.parse(netmhc_software_meta)
    def entry = netmhc_software_meta_map[tool_name]
    // Take OS into account. NetMHC provides different binaries for Mac and Linux
    if (params["netmhc_system"] == 'darwin') {
        entry = netmhc_software_meta_map["${tool_name}_darwin"]
    }
    // If so, add the tool name and user installation path to the external tools import channel
    def ch_netmhc_exe = channel.empty()
    ch_netmhc_exe.bind([
        tool_name,
        entry.version,
        entry.software_md5,
        file(params["${tool_name}_path"], checkIfExists:true),
        entry.data_url ? file(entry.data_url, checkIfExists:true) : [],
        entry.data_md5 ? entry.data_md5 : "",
        entry.binary_name
    ])
    return ch_netmhc_exe
}

// Groovy function to parse JSON and return a list/map
def parseJson(jsonPath) {
    new groovy.json.JsonSlurper().parse(file(jsonPath, checkIfExists: true))
}
