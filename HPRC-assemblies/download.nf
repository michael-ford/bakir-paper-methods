nextflow.enable.dsl=2

process fetchSampleList {

    publishDir "./", mode: 'rellink'

    output:
    path "sample_list.out"

    script:
    """
    curl -L -o sample_list.out 'https://raw.githubusercontent.com/mourisl/T1K_manuscript_evaluation/083d262/HPRC_process/sample_list.out'
    """
}


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
    sample_list_path = fetchSampleList()

    downloadFiles(sample_list_path
                .splitText().map{it -> it.trim()})
}
