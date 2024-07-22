process TCOFFEE_SEQREFORMAT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.46.0.919e8c6b--hfc96bf3_0':
        'biocontainers/t-coffee:13.46.0.919e8c6b--hfc96bf3_0' }"
    storeDir "${params.output}/${aligner}_percid"


    input:
    tuple val(fam), val(aligner), path(infile)

    output:
    tuple path("${prefix}.txt"), emit: formatted_file
    path "versions.yml" , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${fam}"
    """
    export TEMP='./'
    t_coffee -other_pg seq_reformat \
        -in ${infile} \
        -output=sim_idscore \
        > "${prefix}.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${fam}"
    """
    touch "${prefix}.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """
}