#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads    = null
params.outdir   = "results"
params.adapters = null

//----------------------------------------
// Canal: emparejar archivos R1/R2 por nombre de muestra
//----------------------------------------
Channel
    .fromPath(params.reads)
    .ifEmpty { error "No se encontraron FASTQ con el patrón: ${params.reads}" }
    .map { file ->
        def m = file.name =~ /^(.+)_R([12])_001\.fastq\.gz$/
        if (!m.matches()) {
            error "Archivo con nombre inesperado (no encaja con *_R[12]_001.fastq.gz): ${file.name}"
        }
        tuple(m[0][1], m[0][2], file)
    }
    .groupTuple(by: [0])
    .map { sample, reads, files ->
        if (files.size() != 2) {
            error "La muestra ${sample} no tiene exactamente 2 archivos (R1+R2); tiene ${files.size()}: ${files*.name}"
        }
        def sorted = [reads, files].transpose().sort { it[0] }
        tuple(sample, sorted.collect { it[1] })
    }
    .set { fastq_pairs }

fastq_pairs.view()

//----------------------------------------
// Proceso 1: FastQC (control de calidad previo al recorte)
//----------------------------------------
process FASTQC_RAW {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    conda '/home/idalia/miniconda3/envs/rembe_qc'

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
    publishDir "${params.outdir}/trimmed", mode: 'symlink'   // do not copy
    conda '/home/idalia/miniconda3/envs/rembe_qc'

    input:
    tuple val(sample), path(reads)
    path adapters

    output:
    tuple val(sample),
          path("${sample}_R1_paired.fastq.gz"),
          path("${sample}_R2_paired.fastq.gz"), emit: paired
    tuple path("${sample}_R1_unpaired.fastq.gz"),
          path("${sample}_R2_unpaired.fastq.gz"), emit: unpaired

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
    conda '/home/idalia/miniconda3/envs/rembe_qc'

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
    conda '/home/idalia/miniconda3/envs/rembe_qc'

    input:
    path raw_zips,  stageAs: 'fastqc_raw/*'
    path trim_zips, stageAs: 'fastqc_trimmed/*'

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data

    script:
    """
    multiqc fastqc_raw fastqc_trimmed -o .
    """
}

//----------------------------------------
// Workflow principal
//----------------------------------------
workflow {
    // Canal del archivo de adaptadores
    adapters_ch = Channel.fromPath(params.adapters, checkIfExists: true).first()

    // QC previo sobre los RAW
    FASTQC_RAW(fastq_pairs)

    // Recorte
    TRIM(fastq_pairs, adapters_ch)

    // QC posterior sobre las reads recortadas (solo paired)
    FASTQC_TRIMMED(TRIM.out.paired)

    // MultiQC: fusionar informes de QC previo y posterior
    MULTIQC(
        FASTQC_RAW.out.zip.collect(),
        FASTQC_TRIMMED.out.zip.collect()
    )
}