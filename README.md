# Firecloud/WDL DNA methylation workflows
This platform contains publicly accessible cloud-based preprocessing and quality control pipelines that go from raw data to CpG-level methylation estimates. The technologies covered include whole genome bisulfite sequencing (WGBS), reduced representation bisulfite sequencing (RRBS), hybrid selection (capture) bisulfite sequencing (HSBS) and Illumina methylation microarrays. Leveraging the Firecloud platform allows users to: 

1) ensure cross-platform reproducibility of analyses 
2) achieve scalability to large whole genome datasets with 100GB+ of raw data per sample, and to single-cell datasets with thousands of cells 
3) provide access to best-practice analysis pipelines  
4) enable integration and comparison between user-provided data and publicly available data (e.g. TCGA)


### To run the analysis in the FireCloud platform
Analysis should be run in two successive processes: 1) Alignment and methylation calling
													2) Aggregation and quality control analysis

## Alignment and methylation calling
In order to perform alignment and methylation calling choose *bsseq_preprocess_pe_dynamic_disk* method configuration


1) Upload the .fastq files to the google cloud bucket
2) 

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

Note: On Mac, cromwell can be installed with the following command
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
