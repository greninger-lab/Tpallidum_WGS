# Tpallidum WGS
This pipeline is intended for assembly and annotation of Treponema pallidum whole genomes.

This pipeline takes gzipped fastq files and outputs consensus fastas annotated with Prokka. Running on the cloud is recommended due to memory-intensive mapping steps. 

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
- Example command for fastqs in current directory: ```nextflow run greninger-lab/Tpallidum_WGS --INPUT ./ --OUTDIR output/ -r gates_genomes -resume -with-trace -c ~/nextflow.covid.config -profile Cloud```


| Command  | Description |
| ---      | ---         | 
| --INPUT  | Input folder where gzipped fastqs are located. For current  directory, `./` can be used.
| --OUTDIR | Output folder where .bams and consensus fastas will be piped into.
| -resume  | nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-docker ubuntu:18.04 | Runs command with Ubuntu docker.
| -with-trace | Outputs a trace.txt that shows which processes end up in which work/ folders. 
