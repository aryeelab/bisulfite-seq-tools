# Firecloud/WDL DNA methylation workflows
This platform contains publicly accessible cloud-based preprocessing and quality control pipelines that go from raw data to CpG-level methylation estimates. The technologies covered include whole genome bisulfite sequencing (WGBS), reduced representation bisulfite sequencing (RRBS), hybrid selection (capture) bisulfite sequencing (HSBS) and Illumina methylation microarrays. Leveraging the Firecloud platform allows users to: 

1) ensure cross-platform reproducibility of analyses 
2) achieve scalability to large whole genome datasets with 100GB+ of raw data per sample, and to single-cell datasets with thousands of cells 
3) provide access to best-practice analysis pipelines  
4) enable integration and comparison between user-provided data and publicly available data (e.g. TCGA)


## To run the analysis in the FireCloud platform
Analysis should be run in two successive processes: 
1) Alignment and methylation calling
2) Aggregation and quality control analysis

Before running the processes, you need to generate participants file and participant_set file. Both of these files are tab separated text files. Examples of these files are shown in *Firecloud_imports* subdirectory.


## Alignment and methylation calling
In order to perform alignment and methylation calling choose *bismark_rrbs* or *bismark_wgbs* method configuration. As the name indicates
*bismark_rrbs* is for samples that are generated from Reduced Representation Bisulfite Sequencing (RRBS) with Mspl digestion and *bismark_wgbs*
is for data generated from Whole Genome Bisulfite Sequencing (WGBS).

1) Upload the .fastq files to the google cloud bucket
2) In the FireCloud workspace choose *bismark_rrbs* or *bismark_wgbs* method configuration
3) In order to change the reference genome index, click **Edit Configuration** and change the genome_index in the list of inputs
4) Change other parameters according to preference
4) Press *Launch Analysis* in upper right hand corner
5) Choose the participants from the list of files
6) Click **Launch**

You can observe the status of the job by going to *Monitor* tab

## Aggregation and Quality Control Analysis
After the alignment and methylation calling each sample will have their methylation information and metadata stored in HDF5 format
In order to aggregate all of them and obtain the quality control report
1) Choose *aggregate_bismark_output* method configuration
2) Choose the right BSGenome package in the *BSGenome_package* option and choose the location of the tar.gz file in *BSGenome_tagz* option
3) Save it and press *Launch Analysis* 
4) Since this is the aggregation step the entity root type will be participant_set, so you will choose participant_set with participants of your interest
5) Finally click **Launch**

To check the results from any of the work flows, go to Monitor tab, click *View* in the Status columns and then click the *Workflow ID* in the bottom of the page.


## Running the WDL workflow in local environment
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


## FAQ
