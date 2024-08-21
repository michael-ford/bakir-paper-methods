#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${workflow.projectDir}/HPRC-assemblies"  // Base directory for data
params.script_path = "${workflow.projectDir}/Immuannot/scripts/immuannot.sh"
params.resources_dir = "${workflow.projectDir}/Immuannot-data/Data-2023Oct27"
params.publish_directory = "${workflow.projectDir}/HPRC-Skirt-annotations"  // Directory to publish results
params.threads = 12

// Define a process to run immuannot.sh
process RunImmuannot {
    tag "${sample_id}"
    publishDir "${params.publish_directory}", mode: 'rellink'

    input:
        tuple val(sample_id), path(contig_file)

    output:
        path "${sample_id}"

    script:
        """
        bash ${params.script_path} --contig ${contig_file} -r ${params.resources_dir} -o ${sample_id} -t ${params.threads}
        mv ${sample_id}.gtf.gz ${sample_id}/
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

    RunImmuannot(contig_files)
}

// Define a workflow that manages the data and process execution
// workflow {
//     // Testing with a specific file
//     test_file = file("${workflow.projectDir}/HPRC-assemblies/HG00438/HG00438.maternal.f1_assembly_v2.fa.gz")
//     sample_id = 'HG00438.maternal'

//     // Call the process for the test file
//     RunImmuannot(tuple(sample_id, test_file))
// }
