//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    // SAMPLESHEET_CHECK ( samplesheet )
    Channel.fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.sample         = row.sample
    meta.alleles        = row.alleles
    meta.anno           = row.anno

    def array = []
    if (!file(row.filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> file does not exist!\n${row.Filename}"
    } else {
        array = [ meta, file(row.filename) ]
    }
    return array
}
