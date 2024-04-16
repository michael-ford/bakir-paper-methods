#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.skirt_dir = '../skirt' // Default relative path
// Convert relative path to absolute path
String skirtAbsolutePath = new File(params.skirt_dir).getAbsolutePath()


process runSkirt {
    tag "${subject}"

    publishDir "./", mode: 'rellink'

    input:
    tuple val(subject), val(sample), path(file)

    output:
    path "${subject}/*"

    script:
    """
    mkdir -p ${subject}
    cd ${subject}
    echo ${file}
    export SKIRT_WD=${skirtAbsolutePath}
    ${skirtAbsolutePath}/scripts/miniskirt.sh ../${file} ${sample}
    """
}

// Workflow Definition
workflow {
    Channel.of(tuple('HG00438', 'HG00438.paternal.f1_assembly_v2', file('../HPRC-assemblies/HG00438/HG00438.paternal.f1_assembly_v2.fa.gz'))) | runSkirt
}
