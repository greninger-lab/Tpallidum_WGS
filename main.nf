#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    T. pallidum whole genome sequencing
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run michellejlin/Tpallidum_WGS \\
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false

// if INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
   params.INPUT = "${params.INPUT}/"
}
// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}

// Reference files
ADAPTERS = file("${baseDir}/All_adapters.fa")
REF_FASTAS = file("${baseDir}/refs/TPA_refseqs.fasta")
REF_FASTAS_TRIM = file("${baseDir}/refs/TPA_refseqs_trim.fasta")
NC_021508 = file("${baseDir}/refs/NC_021508.fasta")
// bowtie2 indexes
NC_021508_1 = file("${baseDir}/refs/NC_021508.1.bt2")
NC_021508_2 = file("${baseDir}/refs/NC_021508.2.bt2")
NC_021508_3 = file("${baseDir}/refs/NC_021508.3.bt2")
NC_021508_4 = file("${baseDir}/refs/NC_021508.4.bt2")
NC_021508_5 = file("${baseDir}/refs/NC_021508.rev.1.bt2")
NC_021508_6 = file("${baseDir}/refs/NC_021508.rev.2.bt2")
// bwa indexes
NC_021508_BWA1 = file("${baseDir}/refs/NC_021508.fasta.amb")
NC_021508_BWA2 = file("${baseDir}/refs/NC_021508.fasta.ann")
NC_021508_BWA3 = file("${baseDir}/refs/NC_021508.fasta.bwt")
NC_021508_BWA4 = file("${baseDir}/refs/NC_021508.fasta.pac")
NC_021508_BWA5 = file("${baseDir}/refs/NC_021508.fasta.sa")

// Scripts
TP_MAKE_SEQ = file("${baseDir}/tp_make_seq.R")
TP_GENERATE_CONSENSUS = file("${baseDir}/tp_generate_consensus.R")


// Read in fastq pairs into input_read_ch
if(params.SINGLE_END == false){ 
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Adapter trim with Trimmomatic
process trimReads { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), file(R1), file(R2) from input_read_ch
        file ADAPTERS
    output: 
        tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") into Trim_out_ch

    script:
    """
    #!/bin/bash

    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

    """
}

// Use bbduk to filter reads that match Tp genomes
process filterTp {
    container "quay.io/thanhleviet/bbtools:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") from Trim_out_ch
        file REF_FASTAS_TRIM
    output:
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") into Trimmed_filtered_reads_ch1
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") into Trimmed_filtered_reads_ch2
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") into Trimmed_filtered_reads_ch3

    script:
    """
    #!/bin/bash
    bbduk.sh in1='${base}.R1.paired.fastq.gz' in2='${base}.R2.paired.fastq.gz' out1='${base}_unmatched_r1.fastq.gz' out2='${base}_unmatched_r2.fastq.gz' outm1='${base}_matched_r1.fastq.gz' outm2='${base}_matched_r2.fastq.gz' ref=${REF_FASTAS_TRIM} k=31 hdist=2 stats='${base}_stats_tp.txt' overwrite=TRUE t=${task.cpus} -Xmx30g

    """
}

process mapReads {
    container "quay.io/biocontainers/bowtie2:2.4.1--py37h8270d21_3"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") from Trimmed_filtered_reads_ch1
        file(NC_021508)
        file(NC_021508_1)
        file(NC_021508_2)
        file(NC_021508_3)
        file(NC_021508_4)
        file(NC_021508_5)
        file(NC_021508_6)
    output:
        tuple val(base),file("${base}.sam") into Aligned_sam_ch

    script:
    """
    #!/bin/bash

    bowtie2 -x NC_021508 -1 '${base}_matched_r1.fastq.gz' -2 '${base}_matched_r2.fastq.gz' -p ${task.cpus} > ${base}.sam
    """
}

process samToBam {
    container "quay.io/biocontainers/samtools:1.6--h9dace67_6"

    input:
        tuple val(base),file("${base}.sam") from Aligned_sam_ch
    output:
        tuple val(base),file("${base}_firstmap.sorted.bam") into Sorted_bam_ch
    
    script:
    """
    #!/bin/bash

    /usr/local/bin/samtools view -bS ${base}.sam > ${base}.bam
    /usr/local/bin/samtools sort -o ${base}_firstmap.sorted.bam ${base}.bam
    """

}

process deNovoAssembly {
    container "quay.io/biocontainers/unicycler:0.4.4--py37h13b99d1_3"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") from Trimmed_filtered_reads_ch2
    output:
        tuple val(base),file("assembly.gfa"),file("assembly.fasta") into Unicycler_ch
    script:
    """
    #!/bin/bash

    /usr/local/bin/unicycler -1 ${base}_matched_r1.fastq.gz -2 ${base}_matched_r2.fastq.gz -o ./ -t ${task.cpus}

    """
}

// Merges assembly and mapping to make consensus sequence
process mergeAssemblyMapping {
    container "quay.io/michellejlin/tpallidum_wgs:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1

    input:
        tuple val(base),file("assembly.gfa"),file("assembly.fasta") from Unicycler_ch
        tuple val(base),file("${base}_firstmap.sorted.bam") from Sorted_bam_ch
        file(NC_021508)
        file(NC_021508_BWA1)
        file(NC_021508_BWA2)
        file(NC_021508_BWA3)
        file(NC_021508_BWA4)
        file(NC_021508_BWA5)
        file(TP_MAKE_SEQ)
    output:
        tuple val(base),file("${base}_consensus.fasta") into Consensus_ch

    script:
    """
    #!/bin/bash
    echo ${base}

    Rscript --vanilla ${TP_MAKE_SEQ} \'${base}\' \'NC_021508\'
    """
}

process remapReads {
    container "quay.io/michellejlin/tpallidum_wgs"

    input:
        tuple val(base),file("${base}_consensus.fasta") from Consensus_ch
        tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") from Trimmed_filtered_reads_ch3
    output:
        tuple val(base),file("${base}_remapped.sorted.bam") into Remapped_bam_ch

    script:
    """
    /usr/local/miniconda/bin/bowtie2-build -q ${base}_consensus.fasta ${base}_aligned_scaffolds_NC_021508
    /usr/local/miniconda/bin/bowtie2 -x ${base}_aligned_scaffolds_NC_021508 -1 ${base}_matched_r1.fastq.gz -2 ${base}_matched_r2.fastq.gz -p ${task.cpus} | /usr/local/miniconda/bin/samtools view -bS - > ${base}_remapped.bam 
    /usr/local/miniconda/bin/samtools sort -o ${base}_remapped.sorted.bam ${base}_remapped.bam
    """
}

// process generateConsensus {
//     container "quay.io/michellejlin/tpallidum_wgs"

//     input:
//         tuple val(base),file("${base}_remapped.sorted.bam") from Remapped_bam_ch
//         tuple val(base),file("${base}_consensus.fasta") from Consensus_ch
//         tuple val(base),file("${base}_matched_r1.fastq.gz"), file("${base}_matched_r2.fastq.gz") from Trimmed_filtered_reads_ch3
//     output:
//         tuple val(base),file("${base}_remapped.sorted.bam")

//     script:
//     """
//     /usr/local/miniconda/bin/bowtie2-build -q ${base}_consensus.fasta ${base}_aligned_scaffolds_NC_021508
//     /usr/local/miniconda/bin/bowtie2 -x ${base}_aligned_scaffolds_NC_021508 -1 ${base}_matched_r1.fastq.gz -2 ${base}_matched_r2.fastq.gz -p ${task.cpus} | /usr/local/miniconda/bin/samtools view -bS - > ${base}_remapped.bam 
//     /usr/local/miniconda/bin/samtools sort -o ${base}_remapped.sorted.bam ${base}_remapped.bam
//     """

// }


} 

// Single end
else { 
    input_read_ch = Channel
        .fromPath("${params.INPUT}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}

}