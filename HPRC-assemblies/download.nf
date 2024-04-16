nextflow.enable.dsl=2

process downloadFiles {

    publishDir "./", mode: 'rellink'

    input:
    val(line)

    output:
    path("${line}"), emit: sample_dir

    script:
    """
    mkdir -p ${line}
    cd ${line}
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/${line}/assemblies/year1_freeze_assembly_v2/${line}.maternal.f1_assembly_v2.fa.gz
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/${line}/assemblies/year1_freeze_assembly_v2/${line}.paternal.f1_assembly_v2.fa.gz
    """
}

workflow {
    downloadFiles(Channel.fromPath("../T1K-paper-data/Supplemental_Code/T1K_manuscript_evaluation/HPRC_process/sample_list.out")
                .splitText().map{it -> it.trim()})
}
