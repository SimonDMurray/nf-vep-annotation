# code-challenge-nextflow-vep-annotation
A 48 hour challenge to implement vep annotation using Nextflow

## Pre-requisites
1. Ensure you have a local installation of Nextflow, ideally > `22.10.X`
2. Ensure you have a local instllation of Docker

## Building the Docker Image
From the directory that contains the `Dockerfile`, run `docker build -t vep:v0.1 .`

If you choose a to name your image differently, then change the line `container = 'vep:v0.1'` in `nextflow.config` to your image name

## Running the pipeline

There is a copy of the input data at `example-data/subset-sorted.vcf.gz`

Assuming your working directory is this repository and you have a Nextflow binary in you `PATH` the command to run the pipeline should be `nextflow run main.nf --input example-data/subset-sorted.vcf.gz`

## Things I would have liked to add
1. Currently the pipeline only handles `vcf.gz` data, I would like to have included a parameter switch that allowed it to handle `.vcf` or `.vcf.gz`
2. The container is a local image provided, ideally this image would be pushed to DockerHub and then the image would be pulled if there is no local copy
3. The vcf produced is not sorted due to the nature of how Nextflow parallelises tasks, adding `picard` to the container and a process that sorts the output vcf would fix this 
4. The regenerating of the plugin header lines is tool specific and does not scale well as you have to know how each tool updates the header
5. Adding the scratch directive if running on an executor other than local will mean all the little vcf files will be stored on execution node so the file sysystem doesn't get overloaded
