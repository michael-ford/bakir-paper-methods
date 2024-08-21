nextflow.enable.dsl=2

params.data_dir = "${workflow.projectDir}/HPRC-assemblies"  // Base directory for data
params.publish_directory = "${workflow.projectDir}/HPRC-assemblies-annotations"  // Directory to publish results

process runBakir {
    /*
    Process to emulate the Bash script functionality.
    Assumes 'bakir' is available and configured correctly in the environment.
    */
    label 'bakir_process'
    publishDir "${params.publish_directory}", mode: 'rellink'

    input:
    tuple val(subject), path(input_file)

    output:
    path "$subject/*"  // Output everything in the subject directory

    script:
    """
    # Create a directory based on the subject name
    mkdir -p $subject

    # Extract the base name of the file without the extension
    base_name=\$(basename $input_file .fa.gz)

    # Run kir-annotator and redirect output and error
    bakir -o "$subject/\$base_name" "$input_file"
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

    runBakir(contig_files)
}

// workflow {
//     // Collect files based on the glob pattern
//     input_files = Channel.fromPath("./HPRC-assemblies/HG00438/HG00438.maternal.f1_assembly_v2.fa.gz")

//     // Create a channel of tuples containing the subject and the file
//     subject_files = input_files.map { file ->
//         def subject = file.baseName.split(/\.f1_assembly/)[0]  // Adjust regex as needed to extract the subject name
//         return tuple(subject, file)
//     }

//     // Call the main annotation process for each tuple
//     runBakir(subject_files)
// }
