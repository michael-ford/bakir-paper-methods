#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info('Running SKIRT Workflow...')

// Define parameters
params.data_dir = "${workflow.projectDir}/HPRC-assemblies"  // Base directory for data
params.skirt_dir = "${workflow.projectDir}/skirt" // Default relative path
// Convert relative path to absolute path
String skirtAbsolutePath = new File(params.skirt_dir).getAbsolutePath()


process runSkirt {
    tag "${subject}"

    publishDir "${workflow.projectDir}/HPRC-Skirt-annotations", mode: 'rellink'

    input:
    tuple val(sample), path(file)

    output:
    path "${sample}"

    script:
    """
    mkdir -p ${sample}
    export SKIRT_WD=${params.skirt_dir}
    ${params.skirt_dir}/scripts/miniskirt.sh ${file} ${sample}
    """
}

workflow {
    contig_files = Channel.fromFilePairs("${params.data_dir}/*/*.{maternal,paternal}.f1_assembly_v2.fa.gz",
        size: 1,
        flat: true,
        checkIfExists: true
    ).map { sample_id, file ->
        def parts = file.baseName.split("\\.")
        def id = parts.take(parts.size() - 2).join(".")
        return tuple(id, file)
    }

    runSkirt(contig_files)
}


// workflow {
//     // Testing with a specific file
//     test_file = file("${params.data_dir}/HG00438/HG00438.maternal.f1_assembly_v2.fa.gz")
//     def parts = test_file.baseName.split("\\.")
//     def id = parts.take(parts.size() - 2).join(".")

//     // Call the process for the test file
//     runSkirt(tuple(id, test_file))
// }
