nextflow.enable.dsl = 2


// Define the relative path in params
params.skirt_dir = "../skirt"

// Dynamically evaluate the absolute path
String skirtAbsolutePath = new File(params.skirt_dir).getAbsolutePath()

// Import the runSkirt process and pass the absolute path
include { runSkirt } from '../HPRC-Skirt-annotations.nf' params(skirt_dir: skirtAbsolutePath, publish_directory: "./skirt-annotations")

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
    runSkirt(bakir_inputs)
}
