//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true )
        .map { get_samplesheet_paths(it) }
        .set { meta }

    emit: meta                  // channel: [ val(meta), [ files ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, filenames ]
def get_samplesheet_paths(LinkedHashMap row) {
    //---------
    // Save sample, alleles, mhc_class and file_type in a dictionary (metadata)
    // and return a list of meta and the filename.
    //---------

    def meta = [:]
    meta.sample         = row.sample
    meta.alleles        = row.alleles
    meta.mhcclass       = row.mhc_class
    meta.inputtype      = row.inputtype
    expression = row.expression ? file(row.expression, checkIfExists: true) : []

    def array = []
    if (!file(row.filename).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> file does not exist!\n${row.Filename}"
    }
    else {
        array = [meta, file(row.filename)]
    }
    return array
}

def generate_allele_string(String alleles, String mhcclass) {
    // Collect the allele information from the file
    def allele_string
    valid_class1_loci = ['A*','B*','C*','E*','G*']
    valid_class2_loci = ['DR','DP','DQ']
    if ( alleles.endsWith(".txt") || alleles.endsWith(".alleles") )  {
        allele_string = file(alleles).readLines().join(';')
        if ((mhcclass == 'I' & valid_class2_loci.any { allele_string.contains(it)}) |
        (mhcclass == 'II' & valid_class1_loci.any { allele_string.contains(it)})) {
            exit 1, "ERROR: Please check input samplesheet -> invalid mhc class and allele combination found!\n${row.Filename}"
        }
    }
    // or assign the information to a new variable
    else {
        allele_string = alleles
    }
    return allele_string
}
