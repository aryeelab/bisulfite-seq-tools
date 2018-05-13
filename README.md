# Firecloud/WDL DNA methylation workflows
This platform contains publicly accessible cloud-based preprocessing and quality control pipelines that go from raw data to CpG-level methylation estimates. The technologies covered include Whole Genome Bisulfite Sequencing (WGBS), Reduced Representation Bisulfite Sequencing (RRBS) and Hybrid Selection (capture) Bisulfite Sequencing (HSBS). Leveraging the Firecloud platform allows users to:

1) ensure cross-platform reproducibility of analyses
2) achieve scalability to large whole genome datasets with 100GB+ of raw data per sample, and to single-cell datasets with thousands of cells 
3) provide access to best-practice analysis pipelines  
4) enable integration and comparison between user-provided data and publicly available data (e.g. TCGA)


## Workflow steps
Analysis should be run in two successive processes: 
1) Alignment and methylation calling
2) Aggregation and quality control analysis


Before running the processes, you need to generate participants file and participant_set file. Both of these files are tab separated text files. Examples of these files are shown in *Firecloud_imports* subdirectory.


## Alignment and methylation calling
In order to perform alignment and methylation calling choose *bismark_rrbs*, *bismark_wgbs* or *bismark_hsbs* method configuration with appropriate reference genome suffix. As the name indicates
*bismark_rrbs* is for samples that are generated from Reduced Representation Bisulfite Sequencing (RRBS) with Mspl digestion and *bismark_wgbs* is for data generated from Whole Genome Bisulfite Sequencing (WGBS). *bismark_hsbs* is for data generated from Hybrid Selection Bisulfite Sequencing (HSBS). These worflows can also combine fastq files from multiple lanes if the samples are sequenced in such a way.


1) Upload the fastq files to the Google cloud bucket
2) Upload additional files such as target coverage bed file for HSBS sequencing
3) In the FireCloud workspace choose *bismark_rrbs*, *bismark_wgbs* or *bismark_hsbs* method configuration with appropriate reference genome suffix
4) Change other parameters according to preference
5) Press *Launch Analysis* in upper right hand corner
6) Choose the participants from the list of files
7) Click **Launch**


If you are interested in conducting alignment and methylation calling for an entire participant set. Choose the participant set and in the box named **Define expression** type *this.participants*


You can observe the status of the job by going to *Monitor* tab

### Alignment and methylation calling specific parameters
```
r1_fastq & r2_fastq
```
- Paired end FASTQ files

```
samplename
```
- base name for a sample

```
genome_index
```
- Reference genome to conduct bisulfite specicfic alignment

```
n_bp_trim_read1 & n_bp_trim_read1
```
- Number of bases to trim in read 1 and read 2. This is only specified for wgbs and hsbs workflows. rrbs workflow uses the --rrbs option from trimGalore.

```
chrom_sizes
```
- chrom_sizes file in order to generate the BIGWIG file



## Aggregation and Quality Control Analysis
After the alignment and methylation calling each sample will have their methylation information and metadata stored in HDF5 format

In order to aggregate all of them and obtain the quality control report
1) Choose *aggregate_bismark_output* method configuration with appropriate reference genome suffix
2) Change other parameters according to preference
3) Save it and press *Launch Analysis* 
4) Since root entity type in aggregation step is participant set, you will choose participant_set with participants of your interest
5) Finally click **Launch**

To check the results from any of the work flows, go to Monitor tab, click *View* in the Status columns and then click the *Workflow ID* in the bottom of the page.

### Aggregation specific parameters
```
in_pe_reports_files
```
- Report files from the bismark alignment

```
in_covgz_files
```
- coverage output file

```
in_mbias_files
```
- M-bias file

```
BSGenome_targz
```
- BSGenome targz file

```
BSGenome_package
```
- Name of the BSGenome package

```
Genome_build
```
- Name of the genome build, e.g. mm10, hg19 etc. 


### General parameters
```
cpu
```
- Number of CPUs

```
disks
```
- Name of the disk (Depends on the disk size)

```
memory
```
- Memory requirement

```
multicore
```
- Number of cores to use, especially during alignment

```
preemptible
```
- Preemptible option is to use preemptible virtual machine, if it is set to 0 the workflow runs uninterrupted. If it is set to any integer other than zero it may be interpreted that many times. However preemtible option reduces the computing cost significantly.




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
wget https://storage.googleapis.com/aryeelab/bismark-index/mm10_chr19/bismark_mm10_chr19.tar.gz
```

### Running the WDL workflow in Cromwell

Note: On Mac, cromwell can be installed with the following command
```
brew install cromwell
```

Following commands are based on cromwell version 30. Edit the file paths in the json file according to your directory paths to test the workflows.

#### WGBS, Paired-end reads
```
java -jar cromwell-30.2.jar run call_bismark_wgbs.wdl -i sample1_wgbs_pe.json
```

#### RRBS, Paired-end reads
```
java -jar cromwell-30.2.jar run call_bismark_rrbs.wdl -i sample1_rrbs_pe.json
```


#### HSBS, paired-end reads
```
java -jar cromwell-30.2.jar run call_bismark_hsbs.wdl -i sample1_hsbs_pe.json
```

#### Single-end reads (Out of date)
```
java -jar cromwell-30.2.jar run bsseq_preprocess_se.wdl -i sample1_se.json
```


## Questions and Comments
Please use the Github Issue tracker with any issue you face with the platform. Any specific questions or comments, contact me at divyswar01@g.harvard.edu
