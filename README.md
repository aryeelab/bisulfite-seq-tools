# Firecloud/WDL DNA methylation workflows

### Setup

1. Clone this repository

2. Build the docker image

```
cd Docker/bismark
docker build -t aryeelab/bismark .
```

3. If you want to run the examples below, download this small genome index:    
```
gsutil cp gs://fc-dceaadae-be69-41ab-a230-0b735c0556c1/bismark_index/bismark_mm10_chr19.tar.gz testdata/
```

### Running the WDL workflow in Cromwell

Note: On Macs cromwell can be installed with 
```
brew install cromwell
```

#### Paired-end reads
cromwell run bsseq_preprocess_pe.wdl sample1_pe.json

#### Single-end reads
cromwell run bsseq_preprocess_se.wdl sample1_se.json

