//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {

    // Collect the allele information from the file
    def alleleString
    if ( row.alleles.endsWith("txt") )  {
        alleleString = file(row.alleles).readLines().join(';')
    // or assign the information to a new variable
    } else {
        alleleString = row.alleles
    }

    def meta = [:]
    meta.sample         = row.sample
    meta.alleles        = alleleString
    meta.anno           = row.anno
    meta.ext            = row.ext

    def array = []
    if (!file(row.filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> file does not exist!\n${row.Filename}"
    } else {
        array = [ meta, file(row.filename) ]
    }
    return array
}
