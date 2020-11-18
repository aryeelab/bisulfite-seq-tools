workflow bsseq {

  String version = "v1.4"
  
  # 'dev' pipeline versions use the image with the 'latest' tag.
  # release pipeline versions use images tagged with the version number itself
  String image_id = sub(version, "dev", "latest")

  String r1_fastq
  String r2_fastq
  File genome_index
  File monitoring_script
  File chrom_sizes
  String samplename
  String assay_type
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

  scatter (fastq_pair in fastq_pairs) {  
        call trim {input:   r1_fastq = fastq_pair.left, 
                            r2_fastq = fastq_pair.right, 
                            assay_type = assay_type,
                            n_bp_trim_read1 = n_bp_trim_read1, 
                            n_bp_trim_read2 = n_bp_trim_read2, 
                            monitoring_script = monitoring_script, memory = memory, disks = disks, cpu = cpu, preemptible = preemptible, image_id = image_id
        
        }
        call align {input:  r1_fastq = trim.trimmed_r1_fastq, 
                            r2_fastq = trim.trimmed_r2_fastq, 
                            genome_index = genome_index, 
                            samplename = samplename, 
                            multicore = multicore, 
                            monitoring_script = monitoring_script, memory = memory, disks = disks, cpu = cpu, preemptible = preemptible, image_id = image_id
        } 
  }
  
  call merge_replicates {input: bams = align.bam, 
                                reports = align.output_report, 
                                samplename = samplename, 
                                trim_logs = trim.log,
                                align_logs = align.log,
                                assay_type = assay_type,
                                chrom_sizes = chrom_sizes, 
                                genome_index = genome_index, 
                                multicore = multicore, 
                                monitoring_script = monitoring_script, memory = memory, disks = disks, cpu = cpu, preemptible = preemptible, image_id = image_id
  }
  
  


  
  output {
        File bam = merge_replicates.output_bam
        File bai = merge_replicates.output_bai
        File cov_gz = merge_replicates.output_covgz
        File report = merge_replicates.output_report
        File mbias = merge_replicates.mbias_report
        File bigwig = merge_replicates.output_bigwig
        String pipeline_version = version_info.pipeline_version
        String assay = assay_type
        File log = merge_replicates.log
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



task trim {
  String image_id
  File monitoring_script

  File r1_fastq
  File r2_fastq
  String assay_type
  Int n_bp_trim_read1
  Int n_bp_trim_read2
  
  String memory
  String disks
  Int cpu
  Int preemptible

  command {
    chmod u+x ${monitoring_script}
    ${monitoring_script} > monitoring.log &

    cp ${r1_fastq} r1.fastq.gz
    cp ${r2_fastq} r2.fastq.gz
    
    echo "Mode: ${assay_type}" | tee -a log.txt
    
    if [ "${assay_type}" == "rrbs" ] || [ "${assay_type}" == "RRBS" ]; then
        echo "RRBS trimming"  | tee -a log.txt
        trim_galore --rrbs --paired r1.fastq.gz r2.fastq.gz
        mv r1_val_1.fq.gz r1.fastq.gz 
        mv r2_val_2.fq.gz r2.fastq.gz 
    else
        echo "Not applying RRBS trimming" | tee -a log.txt
    fi
    
	echo "5' Trimming options:" | tee -a log.txt
	echo "  n_bp_trim_read1 = ${n_bp_trim_read1}" | tee -a log.txt
	echo "  n_bp_trim_read2 = ${n_bp_trim_read2}" | tee -a log.txt	
	
    if (( ${n_bp_trim_read1} > 0 && ${n_bp_trim_read2} > 0)); then
        echo "Trimming both read 1 and read 2 from the 5' end" | tee -a log.txt
    	trim_galore --paired --clip_R1 ${n_bp_trim_read1} --clip_R2 ${n_bp_trim_read2} r1.fastq.gz r2.fastq.gz
    elif (( ${n_bp_trim_read1} > 0 && ${n_bp_trim_read2} == 0)); then
        echo "Trimming only read 1 from the 5' end" | tee -a log.txt
        trim_galore --paired --clip_R1 ${n_bp_trim_read1} r1.fastq.gz r2.fastq.gz
    elif (( ${n_bp_trim_read1} == 0 && ${n_bp_trim_read2} > 0)); then
        echo "Trimming only read 2 from the 5' end" | tee -a log.txt
        trim_galore --paired --clip_R2 ${n_bp_trim_read2} r1.fastq.gz r2.fastq.gz
    else 
        echo "  Not 5' trimming reads" | tee -a log.txt
        cp r1.fastq.gz r1_val_1.fq.gz
        cp r2.fastq.gz r2_val_2.fq.gz
    fi 
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
    File trimmed_r1_fastq = "r1_val_1.fq.gz"
    File trimmed_r2_fastq = "r2_val_2.fq.gz"    
    File monitoring_log = "monitoring.log"
    File log = "log.txt"
  }

}


task align {

  String image_id
  File r1_fastq
  File r2_fastq
  File genome_index
  File monitoring_script
  String samplename
  
  String multicore
  
  String memory
  String disks
  Int cpu
  Int preemptible

  
  command {
    chmod u+x ${monitoring_script}
    ${monitoring_script} > monitoring.log &
    
    mkdir bismark_index
    echo "Using genome index: `basename ${genome_index}`" | tee -a log.txt
    tar zxvf ${genome_index} -C bismark_index
    echo "Genome index MD5sum: `md5sum ${genome_index}`"  | tee -a log.txt
    
    ln -s ${r1_fastq} r1.fastq.gz
    ln -s ${r2_fastq} r2.fastq.gz

    echo "Starting bismark"            	
	bismark --genome bismark_index --multicore ${multicore} -1 r1.fastq.gz -2 r2.fastq.gz  | tee -a log.txt

    mv *bismark_bt2_pe.bam ${samplename}.bam
    mv *bismark_bt2_PE_report.txt ${samplename}_report.txt
    
    echo "Starting bismark2report"    
    bismark2report --alignment_report ${samplename}_report.txt --output ${samplename}_bismark_report.html  | tee -a log.txt
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
    File log = "log.txt"
  }
  
  
}

task merge_replicates {
  String image_id
  String samplename
  String assay_type
  Array[File] bams
  Array[File] reports
  Array[File] trim_logs
  Array[File] align_logs
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

    echo "### TRIM LOGS ###" | tee -a log.txt
    cat ${sep=' ' trim_logs} | tee -a log.txt
    echo -ne "\n\n### ALIGN LOGS ###\n" | tee -a log.txt
    cat ${sep=' ' align_logs} | tee -a log.txt

    echo -ne "\n\n### MERGE LOGS ###\n" | tee -a log.txt
    
    echo "Using genome index: `basename ${genome_index}`" | tee -a log.txt
    mkdir bismark_index
    tar zxf ${genome_index} -C bismark_index
    
    echo "Running samtools merge" | tee -a log.txt
    samtools merge -n ${samplename}.bam ${sep=' ' bams}
    
    if [ "${assay_type}" == "wgbs" ] || [ "${assay_type}" == "WGBS" ]; then
        echo "Deduplicating..." | tee -a log.txt
        samtools sort -n -o ${samplename}.sorted_by_readname.bam ${samplename}.bam 
        deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
        rm ${samplename}.sorted_by_readname.bam ${samplename}.bam
        mv ${samplename}.sorted_by_readname.deduplicated.bam ${samplename}.bam
    else
        echo "Not deduplicating" | tee -a log.txt
    fi
    
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
            File log = "log.txt"
  }

}
