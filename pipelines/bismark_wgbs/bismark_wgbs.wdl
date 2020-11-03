workflow call_bismark_pool {

  String version = "1.3"
  
  # 'dev' pipeline versions use the image with the 'latest' tag.
  # release pipeline versions use images tagged with the version number itself
  String image_id = sub(version, "dev", "latest")

  String r1_fastq
  String r2_fastq
  File genome_index
  File monitoring_script
  File chrom_sizes
  String samplename
  Int n_bp_trim_read1 = 0
  Int n_bp_trim_read2 = 0
  
  String multicore
  
  String memory
  String disks
  Int cpu
  Int preemptible
  
  call version_info { input: image_id = image_id }
  
  # Split the comma-separated string of fastq file names into arrays
  call split_string_into_array as fastq1 {input: str = r1_fastq}
  call split_string_into_array as fastq2 {input: str = r2_fastq}
  
  Array[Pair[File, File]] fastq_pairs = zip(fastq1.out, fastq2.out)
  
  scatter (fastq_pair in fastq_pairs) { call align_replicates {input: r1_fastq = fastq_pair.left, r2_fastq = fastq_pair.right, n_bp_trim_read1 = n_bp_trim_read1, n_bp_trim_read2 = n_bp_trim_read2, genome_index = genome_index, samplename = samplename, multicore = multicore, monitoring_script = monitoring_script,  memory = memory, disks = disks, cpu = cpu, preemptible = preemptible, image_id = image_id} }

  call merge_replicates {input: bams = align_replicates.bam, reports = align_replicates.output_report, samplename = samplename, chrom_sizes = chrom_sizes, genome_index = genome_index, multicore = multicore, monitoring_script = monitoring_script, memory = memory, disks = disks, cpu = cpu, preemptible = preemptible, image_id = image_id}


  
  output {
        File bam = merge_replicates.output_bam
        File bai = merge_replicates.output_bai
        File cov_gz = merge_replicates.output_covgz
        File report = merge_replicates.output_report
        File mbias = merge_replicates.mbias_report
        File bigwig = merge_replicates.output_bigwig
        String pipeline_version = version_info.pipeline_version
  }
  
}

task version_info {
	String image_id
	command <<<
		cat /VERSION
	>>>
	runtime {
            continueOnReturnCode: false
            docker: "gcr.io/aryeelab/bismark:${image_id}"
            cpu: 1
            memory: "1GB"
        }
	output {
	    String pipeline_version = read_string(stdout())
        }
}


task split_string_into_array {
    String str 
    String arr = "{ARR[@]}"
    command {
        IFS=',' read -ra ARR <<< "${str}"
        for i in "$${arr}"; do
            echo "$i" | tr -d " " >> out
        done
    }

    runtime {
        docker: "debian:stretch"
    }
    
    output {
        Array[String] out = read_lines("out")
    }
}

task align_replicates{

  String image_id
  File r1_fastq
  File r2_fastq
  File genome_index
  File monitoring_script
  String samplename
  Int n_bp_trim_read1
  Int n_bp_trim_read2
  
  String multicore
  
  String memory
  String disks
  Int cpu
  Int preemptible

  
  command {
    chmod u+x ${monitoring_script}
    ${monitoring_script} > monitoring.log &
    WORK_DIR=$PWD
    mkdir tmp
    cd tmp
    mkdir bismark_index
    
    tar zxvf ${genome_index} -C bismark_index
    echo "Genome index MD5sum:"
    md5sum ${genome_index}
    
    ln -s ${r1_fastq} r1.fastq.gz
    ln -s ${r2_fastq} r2.fastq.gz
	
	echo "Trimming options:"
	echo "  n_bp_trim_read1 = ${n_bp_trim_read1}"
	echo "  n_bp_trim_read2 = ${n_bp_trim_read2}"	
	
    if (( ${n_bp_trim_read1} > 0 && ${n_bp_trim_read2} > 0)); then
        echo "Trimming both read 1 and read 2 from the 5' end"
    	trim_galore --paired --clip_R1 ${n_bp_trim_read1} --clip_R2 ${n_bp_trim_read2} r1.fastq.gz r2.fastq.gz
    elif (( ${n_bp_trim_read1} > 0 && ${n_bp_trim_read2} == 0)); then
        echo "Trimming only read 1 from the 5' end"
        trim_galore --paired --clip_R1 ${n_bp_trim_read1} r1.fastq.gz r2.fastq.gz
    elif (( ${n_bp_trim_read1} == 0 && ${n_bp_trim_read2} > 0)); then
        echo "Trimming only read 2 from the 5' end"
        trim_galore --paired --clip_R2 ${n_bp_trim_read2} r1.fastq.gz r2.fastq.gz
    else
        echo "Not trimming reads"
        cp r1.fastq.gz r1_val_1.fq.gz
        cp r2.fastq.gz r2_val_2.fq.gz
    fi
            	
	bismark --genome bismark_index --multicore ${multicore} -1 r1_val_1.fq.gz -2 r2_val_2.fq.gz
    mv *bismark_bt2_pe.bam $WORK_DIR/${samplename}.bam
    mv *bismark_bt2_PE_report.txt $WORK_DIR/${samplename}_report.txt
    
    bismark2report --alignment_report $WORK_DIR/${samplename}_report.txt --output $WORK_DIR/${samplename}_bismark_report.html  
  
  }

  runtime {
  	continueOnReturnCode: false
	docker: "gcr.io/aryeelab/bismark:${image_id}"
	memory: memory
    disks: disks
    cpu: cpu
    preemptible: preemptible
 }
  
  output {
    File bam = "${samplename}.bam"
    File output_report = "${samplename}_report.txt"
    File bismark_report_html = "${samplename}_bismark_report.html"
    File monitoring_log = "monitoring.log"
  }
  
  
}

task merge_replicates {
  String image_id
  String samplename
  Array[File] bams
  Array[File] reports
  String multicore
  File monitoring_script
  File chrom_sizes
  File genome_index

  String memory
  String disks
  Int cpu
  Int preemptible
  
  command {
    
    # The file renaming below is necessary since this version of bismark doesn't allow the 
    # use of --multicore with --basename
    chmod u+x ${monitoring_script}
    ${monitoring_script} > monitoring.log &
    mkdir bismark_index
    tar zxf ${genome_index} -C bismark_index
    
    samtools merge -n ${samplename}.bam ${sep=' ' bams}
    
    samtools sort -n -o ${samplename}.sorted_by_readname.bam ${samplename}.bam 
    deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
    rm ${samplename}.sorted_by_readname.bam ${samplename}.bam
    mv ${samplename}.sorted_by_readname.deduplicated.bam ${samplename}.bam
    
    samtools sort -o ${samplename}.sorted.bam ${samplename}.bam
    samtools index ${samplename}.sorted.bam ${samplename}.sorted.bai
              
    bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${samplename}.bam
    
    gunzip "${samplename}.bedGraph.gz"
    
    bedtools sort -i ${samplename}.bedGraph > ${samplename}.sorted.bedGraph        
    /usr/local/bin/bedGraphToBigWig ${samplename}.sorted.bedGraph ${chrom_sizes} ${samplename}.bw
    
    echo """${sep="\n" reports}""" > all_reports.txt
    Rscript --vanilla /pool_report_files.R -i all_reports.txt -s ${samplename}
  }
  
  runtime {
    continueOnReturnCode: false
    docker: "gcr.io/aryeelab/bismark:${image_id}"
    memory: memory
    disks: disks
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
            File monitoring_log = "monitoring.log" 
            File output_bam = "${samplename}.sorted.bam"
            File output_covgz = "${samplename}.bismark.cov.gz"
            File output_report = "${samplename}_report.txt"
            File mbias_report = "${samplename}.M-bias.txt"
            File output_bigwig = "${samplename}.bw"
            File output_bai = "${samplename}.sorted.bai"
  }

}
<<<<<<< HEAD
=======

>>>>>>> releases/v1.3
