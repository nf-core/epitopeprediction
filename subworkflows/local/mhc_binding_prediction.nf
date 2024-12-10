//
// Check input samplesheet and get read channels
//
include { PREPARE_PREDICTION_INPUT } from '../../modules/local/prepare_prediction_input'
include { SYFPEITHI } from '../../modules/local/syfpeithi'
include { MHCFLURRY } from '../../modules/local/mhcflurry'
include { MHCNUGGETS } from '../../modules/local/mhcnuggets'
include { NETMHCPAN } from '../../modules/local/netmhcpan'
include { NETMHCIIPAN } from '../../modules/local/netmhciipan'
include { parse_netmhc_params } from '../../subworkflows/local/utils_nfcore_epitopeprediction_pipeline'
include { EXTERNAL_TOOLS_IMPORT as NETMHCPAN_IMPORT;
        EXTERNAL_TOOLS_IMPORT as NETMHCIIPAN_IMPORT} from '../../modules/local/external_tools_import'
include { MERGE_PREDICTIONS } from '../../modules/local/merge_predictions'

workflow MHC_BINDING_PREDICTION {
    take:
        ch_peptides

    main:
        ch_versions = Channel.empty()
        ch_binding_predictors_out = Channel.empty()

        //TODO: Add //"pattern": "^(syfpeithi|mhcnuggets|mhcflurry|netmhcpan|netmhciipan)(,(syfpeithi|mhcnuggets|mhcflurry|netmhcpan|netmhciipan)){0,4}$",
        // to nextflow_schema.json once this subworkflow is completed

        //prepare the input file
        PREPARE_PREDICTION_INPUT( ch_peptides )
            .prepared
            .transpose()
            .branch {
                meta, file ->
                    syfpeithi : file.name.contains('syfpeithi')
                        return [meta, file]
                    mhcflurry : file.name.contains('mhcflurry')
                        return [meta, file]
                    mhcnuggets : file.name.contains('mhcnuggets')
                        return [meta, file]
                    netmhcpan: file.name.contains('netmhcpan')
                        return [meta, file]
                    netmhciipan: file.name.contains('netmhciipan')
                        return [meta, file]
                    }
            .set{ ch_prediction_input }

        SYFPEITHI ( ch_prediction_input.syfpeithi )
        ch_versions = ch_versions.mix(SYFPEITHI.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(SYFPEITHI.out.predicted)

        MHCFLURRY ( ch_prediction_input.mhcflurry )
        ch_versions = ch_versions.mix(MHCFLURRY.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCFLURRY.out.predicted)

        MHCNUGGETS ( ch_prediction_input.mhcnuggets )
        ch_versions = ch_versions.mix(MHCNUGGETS.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCNUGGETS.out.predicted)

        if ( "netmhcpan" in params.tools.tokenize(",") )
        {
            //TODO: Refactor to support only one netmhc version
            //TODO: Fix parsing version
            NETMHCPAN_IMPORT( parse_netmhc_params("netmhcpan", "4.1") )
            //TODO: Update netmhc container
            //TODO Fix parsing version
            NETMHCPAN ( ch_prediction_input.netmhcpan.combine(NETMHCPAN_IMPORT.out.nonfree_tools) )
            ch_versions = ch_versions.mix(NETMHCPAN.out.versions)
            ch_binding_predictors_out = ch_binding_predictors_out.mix(NETMHCPAN.out.predicted)
        }

        if ( "netmhciipan" in params.tools.tokenize(",") )
        {
            // TODO: External tools import for netmhciipan
            NETMHCIIPAN_IMPORT( parse_netmhc_params("netmhciipan", "4.3") )
            //TODO: Update netmhc container
            NETMHCIIPAN ( ch_prediction_input.netmhcpan.combine(NETMHCIIPAN_IMPORT.out.nonfree_tools) )
            ch_versions = ch_versions.mix(NETMHCIIPAN.out.versions)
            ch_binding_predictors_out = ch_binding_predictors_out.mix(NETMHCIIPAN.out.predicted)
        }
    ch_binding_predictors_out.groupTuple().view()
    MERGE_PREDICTIONS (ch_binding_predictors_out.groupTuple())
    ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    emit:
    predicted = MERGE_PREDICTIONS.out.merged
    versions = ch_versions
}
