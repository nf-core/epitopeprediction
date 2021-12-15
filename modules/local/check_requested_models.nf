process CHECK_REQUESTED_MODELS {

    conda (params.enable_conda ? "bioconda::fred2=2.0.7 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0"
    }

    input:
        tuple val(alleles), path(input_file)
        path(software_versions)

    output:
        path '*.txt', emit: txt // model_report.txt
        path '*.log', emit: log // model_warnings.log
        path "versions.yml", emit: versions


    script:
        def argument = task.ext.args

        if (argument.contains("peptides") == true) {
            argument += " ${input_file}"
        }

        if (params.multiqc_title) {
            argument += "--title \"$params.multiqc_title\""
        }

        """
        check_requested_models.py ${argument} \
            --alleles '${alleles.join(';')}' \
            --mhcclass ${params.mhc_class} \
            --versions ${software_versions} > model_warnings.log

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)"))
            fred2: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('Fred2').version)"))
        END_VERSIONS
        """

}
