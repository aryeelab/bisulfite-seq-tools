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

#### WGBS, Paired-end reads
```
cromwell run call_bismark_wgbs.wdl sample1_wgbs_pe.json
```

#### RRBS, Paired-end reads
```
cromwell run call_bismark_rrbs.wdl sample1_rrbs_pe.json
```


#### Hybrid selection (capture) bisulfite sequencing, paired-end reads
```
cromwell run call_bismark_hsbs.wdl sample1_hsbs_pe.json
```

#### Single-end reads (Out of date)
```
cromwell run bsseq_preprocess_se.wdl sample1_se.json
```
