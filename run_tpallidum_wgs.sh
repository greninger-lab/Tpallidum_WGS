#!/bin/bash

set -e

# In this example, I cloned the GitHub repo FredHutch/nextflow-aws-batch-squared to the following path.
# Set this next line to point to the location you cloned that repo
RUN_HEADLESS_PY=/Users/gerbix/Documents/pipelines/clomp_run_scripts/nextflow-aws-batch-squared-master/run.py

# This just refers to the profile in your local ~/.aws/credentials file
AWS_PROFILE=default

# This refers to the AWS Batch setup
JOB_QUEUE=multi-az

####################
# WORKFLOW OPTIONS #
####################

# Unique Name for this particular batch of samples
# (no spaces, no special characters, keep it simple)
SAMPLE_BATCH_NAME=all_cov_250-1000
INPUT="s3://covid19-input-data/tpallidum_wgs/2020-10-21/all_cov_250-1000/"
# Make sure the OUTDIR always ends with a '/'
OUTDIR="s3://covid19-input-data/tpallidum_wgs/2020-10-21/all_cov_250-1000/2020-12-16_extra_filtering_output/"

# GitHub repo to run
WORKFLOW=michellejlin/Tpallidum_WGS
REVISION=master

# The options below will likely not change very often
ADAPTERS="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/All_adapters.fa"
REF_FASTAS="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/TPA_rRNA_refs.fasta"
REF_FASTAS_MASKED="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/TP_refs_rRNA_masked.fasta"
REF_FASTAS_TRIM="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/TPA_refseqs_trim.fasta"
NC_021508="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta"
REPEAT_FILTER_FASTA="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/repeat_filter_UPDATE.fasta"

# bowtie2 indexes
NC_021508_1="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.1.bt2"
NC_021508_2="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.2.bt2"
NC_021508_3="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.3.bt2"
NC_021508_4="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.4.bt2"
NC_021508_5="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.rev.1.bt2"
NC_021508_6="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.rev.2.bt2"

# bwa indexes
NC_021508_BWA1="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta.amb"
NC_021508_BWA2="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta.ann"
NC_021508_BWA3="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta.bwt"
NC_021508_BWA4="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta.pac"
NC_021508_BWA5="s3://clomp-reference-data/tool_specific_data/Tpallidum_WGS/refs/NC_021508.fasta.sa"

CONFIG_FILE=s3://clomp-reference-data/tool_specific_data/config_files/tpallidum_wgs.config
WORKING_DIRECTORY=s3://covid19-work/

######################
# / WORKFLOW OPTIONS #
######################

# Make a params JSON file
cat > ${SAMPLE_BATCH_NAME}.params.json << EOF
{
    "INPUT": "${INPUT}",
    "OUTDIR": "${OUTDIR}",
    "ADAPTERS": "${ADAPTERS}",
    "REF_FASTAS": "${REF_FASTAS}",
    "REF_FASTAS_TRIM": "${REF_FASTAS_TRIM}",
    "NC_021508": "${NC_021508}",
    "NC_021508_1": "${NC_021508_1}",
    "NC_021508_2": "${NC_021508_2}",
    "NC_021508_3": "${NC_021508_3}",
    "NC_021508_4": "${NC_021508_4}",
    "NC_021508_5": "${NC_021508_5}",
    "NC_021508_6": "${NC_021508_6}",
    "NC_021508_BWA1": "${NC_021508_BWA1}",
    "NC_021508_BWA2": "${NC_021508_BWA2}",
    "NC_021508_BWA3": "${NC_021508_BWA3}",
    "NC_021508_BWA4": "${NC_021508_BWA4}",
    "NC_021508_BWA5": "${NC_021508_BWA5}",
    "REPEAT_FILTER_FASTA": "${REPEAT_FILTER_FASTA}"
}
EOF

# Copy the params JSON to S3
AWS_PROFILE=$AWS_PROFILE aws s3 cp ${SAMPLE_BATCH_NAME}.params.json ${OUTDIR}
# Copy the config file to S3
AWS_PROFILE=$AWS_PROFILE aws s3 cp ${CONFIG_FILE} ${OUTDIR}nextflow.config

AWS_PROFILE=$AWS_PROFILE python3 ${RUN_HEADLESS_PY} \
    --workflow ${WORKFLOW} \
    --revision ${REVISION} \
    --working-directory ${WORKING_DIRECTORY} \
    --name HEAD_${SAMPLE_BATCH_NAME} \
    --config-file ${OUTDIR}nextflow.config \
    --params-file ${OUTDIR}${SAMPLE_BATCH_NAME}.params.json \
    --temporary-volume /var/lib/docker \
    --job-queue ${JOB_QUEUE} \
    --restart-uuid ${SAMPLE_BATCH_NAME} \
    --with-report ${OUTDIR}report.html \
    --with-trace ${OUTDIR}trace.txt 
