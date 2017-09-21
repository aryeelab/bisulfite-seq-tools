# Firecloud/WDL DNA methylation workflows

### Setup

1. Clone this repository
2. If you want to run the examples below, download this small genome index:

    gsutil cp gs://fc-dceaadae-be69-41ab-a230-0b735c0556c1/bismark_index/bismark_mm10_chr19.tar.gz testdata/

### Build docker image
cd Docker/bismark
docker build -t aryeelab/bismark .

### Run WDL in Cromwell

Note: On Macs this can be installed with `brew install cromwell`

#### Paired-end reads
cromwell run bsseq_preprocess_pe.wdl sample1_pe.json

#### Single-end reads
cromwell run bsseq_preprocess_pe.wdl sample1_pe.json

