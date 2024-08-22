nextflow.enable.dsl = 2


// Define parameters
// params.gb_zip = "gb.zip"
params.publish_directory = "./genbank-samples"

// Importing and using the runBakir process from another script
include { runBakir } from '../HPRC-assemblies-annotations.nf' params(publish_directory: "./bakir-annotations")

// // Process to extract gb.zip
// process unpack_genbank_samples {
//     publishDir "genbank-samples", mode: 'rellink'

//     input:
//     path gb_zip

//     output:
//     path "genbank-samples/*.fa" 

//     script:
//     """
//     mkdir -p genbank-samples
//     unzip $gb_zip -d genbank-samples
//     """
// }

params.sample_dir = "./sample-sequences"

// Process to extract sample ID from .fa files
process get_subject_name {
    input:
    path fa_file

    output:
    tuple env(subject), path(fa_file)

    script:
    """
    # Extracting the subject/sample ID from the filename
    subject=\$(basename ${fa_file} .fa)
    """
}


// Define the main workflow
workflow {

    // // Unzip the gb.zip file and extract .fa files
    // genbank_fa_files = unpack_genbank_samples(Channel.fromPath(params.gb_zip))

    // Get the subject name and path of each .fa file
    // bakir_inputs = get_subject_name(Channel.fromPath("${params.sample_dir}/*.fa"))
    bakir_inputs = get_subject_name(Channel.fromPath("${params.sample_dir}/*.fa"))

    // Run the runBakir process on each extracted .fa file
    runBakir(bakir_inputs)
}
