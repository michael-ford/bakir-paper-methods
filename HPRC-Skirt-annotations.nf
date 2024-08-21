#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info('Running SKIRT Workflow...')

// Define parameters
params.data_dir = "${workflow.projectDir}/HPRC-assemblies"  // Base directory for data
params.skirt_dir = "${workflow.projectDir}/skirt" // Default relative path
params.publish_directory = "${workflow.projectDir}/HPRC-Skirt-annotations"  // Directory to publish results
// Convert relative path to absolute path
String skirtAbsolutePath = new File(params.skirt_dir).getAbsolutePath()


process runSkirt {
    errorStrategy 'ignore'
    tag "${subject}"

    publishDir "${publish_directory}", mode: 'rellink'

    input:
    tuple val(sample), path(contig_file)

    output:
    path "${sample}"

    script:
    """
    echo "Running SKIRT for ${sample}..."
    mkdir -p ${sample}
    export SKIRT_WD=${params.skirt_dir}
    ${params.skirt_dir}/scripts/miniskirt.sh ${contig_file} ${sample}
    """
}

// process test {
//     tag "${subject}"

//     publishDir "${workflow.projectDir}/HPRC-Skirt-annotations", mode: 'rellink'

//     input:
//     tuple val(sample), path(contig_file)

//     // output:
//     // path "${sample}"

//     script:
//     """
//     echo "TEST"
//     mkdir -p ${sample}
//     export SKIRT_WD=${params.skirt_dir}
//     ${params.skirt_dir}/scripts/miniskirt.sh ${contig_file} ${sample}
//     """
// }


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
    contig_files.view()

    runSkirt(contig_files)
}


// workflow {
//     // Testing with a specific file
//     test_file = file("${params.data_dir}/HG00438/HG00438.maternal.f1_assembly_v2.fa.gz")
//     def parts = test_file.baseName.split("\\.")
//     def id = parts.take(parts.size() - 2).join(".")

//     println("Test file: ${test_file}")

//     // Call the process for the test file
//     // runSkirt(tuple(id, test_file))
//     test(tuple(id, test_file))
// }
