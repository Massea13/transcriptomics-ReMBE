#!/usr/bin/env nextflow

/*
 * Quality assessment de transcriptomas:
 *   raw FASTQ  ->  FastQC  ->  MultiQC
 */

nextflow.enable.dsl=2

//----------------------------------------
// Parámetros
//----------------------------------------
params.reads    = null          // p.ej. "/path/data/*.fastq.gz"
params.outdir   = "results"
params.threads  = 4

//----------------------------------------
// Perfiles: conda (local) y cluster (SLURM)
//----------------------------------------
profiles {
    conda {
        conda.enabled = true
    }
    cluster {
        process.executor = 'slurm'
        process.queue    = 'short'
        process.cpus     = 2
        process.memory   = '4 GB'
        process.time     = '1h'
        conda.enabled    = true
    }
    local {
        process.executor = 'local'
        conda.enabled    = false
    }
}

//----------------------------------------
// Canales
//----------------------------------------
Channel
    .fromPath(params.reads)
    .ifEmpty { error "No se encontraron FASTQ con el patrón: ${params.reads}" }
    .map { file -> tuple(file.baseName.replaceAll(/_R?[12](_001)?$/, ''), file) }
    .groupTuple()
    .set { fastq_pairs }

//----------------------------------------
// Proceso 1: FastQC
//----------------------------------------
process FASTQC {

    tag "${sample}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(sample), path(fastqs)

    output:
    path "*.zip"  , emit: zip
    path "*.html" , emit: html

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${fastqs}
    """
}

//----------------------------------------
// Proceso 2: MultiQC
//----------------------------------------
process MULTIQC {

    tag "all_samples"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    conda 'bioconda::multiqc=1.21'

    input:
    path '*',  stageAs: 'fastqc_reports/*'
    path '*.zip'

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data"       , emit: data

    script:
    """
    multiqc fastqc_reports -o .
    """
}

//----------------------------------------
// Workflow principal
//----------------------------------------
workflow {

    FASTQC(fastq_pairs)

    MULTIQC(
        FASTQC.out.html.collect(),
        FASTQC.out.zip.collect()
    )
}