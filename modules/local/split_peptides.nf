// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_PEPTIDES {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'splitted', meta:[:], publish_by_meta:[]) }

    input:
        tuple val(meta), path(peptide)

    output:
        tuple val(meta), path("*.tsv"), emit: splitted

    script:
        def prefix = options.suffix ? "${peptide.baseName}_${options.suffix}" : "${peptide.baseName}"

        """
        split_peptides.py --input ${peptide} \\
        --output_base "${prefix}" \\
        $options.args
        """
}
