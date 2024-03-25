process evaluate_nirmsd {
    tag"${id}"
    storeDir "${params.main_dir}/${aligner}_NiRMSD"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0':
        'biocontainers/mulled-v2-a76a981c07359a31ff55b9dc13bd3da5ce1909c1:84c8f17f1259b49e2f7783b95b7a89c6f2cb199e-0' }"

    input:
        tuple val(id), val(aligner), file(msa), file (template), file(pdb) 

    output:
        tuple val(id), file("*.NiRMSD"), emit: nirmsd
    
    script:
    """
    t_coffee -other_pg irmsd ${msa} -template_file ${template} > ${id}_${aligner}.NiRMSD
    """

}