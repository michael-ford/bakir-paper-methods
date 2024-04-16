#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info('Running SKIRT Workflow...')

// Define parameters
params.data_dir = '../HPRC-assemblies'
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
    export SKIRT_WD=${skirtAbsolutePath}
    ${skirtAbsolutePath}/scripts/miniskirt.sh ../${file} ${sample}
    """
}

// Workflow Definition
workflow {

    // Input Gathering
    Path dataDir = file(params.data_dir)
    Channel
        .fromFilePairs("${dataDir}/*/*.fa.gz", size: 1, flat: true)
        .set { faGzFiles }

    // Workflow Execution
    faGzFiles
        .map { subject, file -> tuple(subject, file.baseName.replace(".fa.gz", ""), file) }
        .set { skirtInput }

    runSkirt(skirtInput)
}
