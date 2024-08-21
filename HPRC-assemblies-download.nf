nextflow.enable.dsl=2

process fetchSampleList {
    output:
    path "HPRC_sample_list.out", emit: HPRC_sample_list
    path "HPRC_PLUS_sample_list.out", emit: HPRC_PLUS_sample_list

    script:
    """
    curl -L -o sample_list.out 'https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/b75b54a/sample_metadata/hprc_year1_sample_metadata.txt'
    tail -n+2 < sample_list.out| grep "HPRC_PLUS" | cut -f1 > HPRC_PLUS_sample_list.out
    tail -n+2 < sample_list.out| grep -v "HPRC_PLUS"  | cut -f1 > HPRC_sample_list.out

    """
}

process downloadFiles {
    publishDir "${workflow.projectDir}/HPRC-assemblies", mode: 'rellink'
    errorStrategy 'ignore'

    input:
    val line

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
process downloadFilesPlus {
    publishDir "${workflow.projectDir}/HPRC-assemblies", mode: 'rellink'
    errorStrategy 'ignore'

    input:
    val line

    output:
    path("${line}"), emit: sample_dir

    script:
    """
    mkdir -p ${line}
    cd ${line}
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/${line}/assemblies/year1_freeze_assembly_v2/${line}.maternal.f1_assembly_v2.fa.gz
    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/${line}/assemblies/year1_freeze_assembly_v2/${line}.paternal.f1_assembly_v2.fa.gz
    """
}

workflow {
    sample_list_out = fetchSampleList()

    sample_list_out.HPRC_sample_list
        .map { file -> file.splitText() }
        .flatten()
        .map { it.trim() }
        .filter { it }
        .set { HPRC_sample_list }
    HPRC_sample_list.view()
    
    sample_list_out.HPRC_PLUS_sample_list
        .map { file -> file.splitText() }
        .flatten()
        .map { it.trim() }
        .filter { it }
        .set { HPRC_PLUS_sample_list }


    downloadFiles(HPRC_sample_list)
    downloadFilesPlus(HPRC_PLUS_sample_list)
}
