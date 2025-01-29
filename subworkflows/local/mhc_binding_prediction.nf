include { PREPARE_PREDICTION_INPUT } from '../../modules/local/prepare_prediction_input/main'
include { SYFPEITHI } from '../../modules/local/syfpeithi'
include { MHCFLURRY } from '../../modules/local/mhcflurry'
include { MHCNUGGETS;
        MHCNUGGETS as MHCNUGGETSII} from '../../modules/local/mhcnuggets'
include { NETMHCPAN } from '../../modules/local/netmhcpan'
include { NETMHCIIPAN } from '../../modules/local/netmhciipan'
include { json2map; parse_netmhc_params } from '../../subworkflows/local/utils_nfcore_epitopeprediction_pipeline'
include { EXTERNAL_TOOLS_IMPORT as NETMHCPAN_IMPORT;
        EXTERNAL_TOOLS_IMPORT as NETMHCIIPAN_IMPORT} from '../../modules/local/external_tools_import'
include { MERGE_PREDICTIONS } from '../../modules/local/merge_predictions/main'

workflow MHC_BINDING_PREDICTION {
    take:
        ch_peptides
        supported_alleles_json

    main:
        ch_versions = Channel.empty()
        ch_binding_predictors_out = Channel.empty()

        //prepare the input file
        PREPARE_PREDICTION_INPUT( ch_peptides, supported_alleles_json)
            .prepared
            .transpose()
            .branch {
                meta, json, file ->
                    def allele_input_dict = json2map(json)
                    mhcflurry : (file.name.contains('mhcflurry_input') && allele_input_dict['mhcflurry'])
                        return [meta + [alleles_supported: allele_input_dict['mhcflurry']], file]
                    mhcnuggets : (file.name.contains('mhcnuggets_input') && allele_input_dict['mhcnuggets'])
                        return [meta + [alleles_supported: allele_input_dict['mhcnuggets']], file]
                    mhcnuggetsii : (file.name.contains('mhcnuggetsii_input') && allele_input_dict['mhcnuggetsii'])
                        return [meta + [alleles_supported: allele_input_dict['mhcnuggetsii']], file]
                    netmhcpan: (file.name.contains('netmhcpan_input') && allele_input_dict['netmhcpan-4.1'])
                        return [meta + [alleles_supported: allele_input_dict['netmhcpan-4.1']], file]
                    netmhciipan: (file.name.contains('netmhciipan_input') && allele_input_dict['netmhciipan-4.3'])
                        return [meta + [alleles_supported: allele_input_dict['netmhciipan-4.3']], file]
                    }
            .set{ ch_prediction_input }

        MHCFLURRY ( ch_prediction_input.mhcflurry )
        ch_versions = ch_versions.mix(MHCFLURRY.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCFLURRY.out.predicted)

        MHCNUGGETS ( ch_prediction_input.mhcnuggets )
        ch_versions = ch_versions.mix(MHCNUGGETS.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCNUGGETS.out.predicted)

        MHCNUGGETSII ( ch_prediction_input.mhcnuggetsii )
        ch_versions = ch_versions.mix(MHCNUGGETSII.out.versions)
        ch_binding_predictors_out = ch_binding_predictors_out.mix(MHCNUGGETSII.out.predicted)

        if ( "netmhcpan-4.1" in params.tools.tokenize(",") )
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

        if ( "netmhciipan-4.3" in params.tools.tokenize(",") )
        {
            // TODO: External tools import for netmhciipan
            NETMHCIIPAN_IMPORT( parse_netmhc_params("netmhciipan", "4.3") )
            //TODO: Update netmhc container
            NETMHCIIPAN ( ch_prediction_input.netmhciipan.combine(NETMHCIIPAN_IMPORT.out.nonfree_tools) )
            ch_versions = ch_versions.mix(NETMHCIIPAN.out.versions)
            ch_binding_predictors_out = ch_binding_predictors_out.mix(NETMHCIIPAN.out.predicted)
        }

    MERGE_PREDICTIONS( ch_binding_predictors_out.map {meta, file -> [meta.subMap('sample','alleles','mhc_class','input_type'), file] }.groupTuple())
    ch_versions = ch_versions.mix(MERGE_PREDICTIONS.out.versions)

    emit:
    predicted = MERGE_PREDICTIONS.out.merged
    versions = ch_versions
}
