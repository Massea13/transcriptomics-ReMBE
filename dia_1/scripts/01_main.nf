#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads    = null
params.outdir   = "results"
params.adapters = null   // ruta al archivo de adaptadores, por ejemplo: /ruta/adapters.fa

//----------------------------------------
// Canal: emparejar archivos R1/R2 por nombre de muestra
//----------------------------------------
Channel
    .fromFilePairs(params.reads, size: 2, checkIfExists: true)
    .ifEmpty { error "No se encontraron pares FASTQ con el patrón: ${params.reads}" }
    .set { fastq_pairs }

fastq_pairs.view()   // depuración temporal

//----------------------------------------
// Proceso 1: FastQC (control de calidad previo al recorte)
//----------------------------------------
process FASTQC_RAW {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(sample), path(fastqs)

    output:
    path "*.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${fastqs}
    """
}

//----------------------------------------
// Proceso 2: Trimmomatic (recorte)
//----------------------------------------
process TRIM {
    tag "${sample}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    conda 'bioconda::trimmomatic=0.39'

    input:
    tuple val(sample), path(reads)
    path adapters

    output:
    tuple val(sample), path("${sample}_R1_paired.fastq.gz"), path("${sample}_R2_paired.fastq.gz"), emit: paired
    path "${sample}_R*_unpaired.fastq.gz", emit: unpaired

    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${sample}_R1_paired.fastq.gz ${sample}_R1_unpaired.fastq.gz \\
        ${sample}_R2_paired.fastq.gz ${sample}_R2_unpaired.fastq.gz \\
        ILLUMINACLIP:${adapters}:2:30:10 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:36
    """
}

//----------------------------------------
// Proceso 3: FastQC (control de calidad posterior al recorte)
//----------------------------------------
process FASTQC_TRIMMED {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path "*.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}
    """
}

//----------------------------------------
// Proceso 4: MultiQC (resumen final, incluye QC previo y posterior)
//----------------------------------------
process MULTIQC {
    tag "all_samples"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda 'bioconda::multiqc=1.21'

    input:
    path raw_zips
    path trim_zips

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data

    script:
    """
    multiqc . -o .
    """
}

//----------------------------------------
// Workflow principal
//----------------------------------------
workflow {
    // Canal del archivo de adaptadores
    adapters_ch = Channel.fromPath(params.adapters, checkIfExists: true)

    // QC previo
    FASTQC_RAW(fastq_pairs)

    // Recorte
    TRIM(fastq_pairs, adapters_ch)

    // QC posterior
    FASTQC_TRIMMED(TRIM.out.paired)

    // MultiQC: fusionar informes de QC previo y posterior
    MULTIQC(
        FASTQC_RAW.out.zip.collect(),
        FASTQC_TRIMMED.out.zip.collect()
    )
}