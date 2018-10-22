workflow call_bismark_pool {

  String r1_fastq
  String r2_fastq
  File genome_index
  File monitoring_script
  File chrom_sizes
  String samplename
  Int n_bp_trim_read1
	Int n_bp_trim_read2
  
  String multicore
  
  String memory
  String disks
  Int cpu
  Int preemptible
  
  
  
   # Split the comma-separated string of fastq file names into arrays
  call split_string_into_array as fastq1 {input: str = r1_fastq}
  call split_string_into_array as fastq2 {input: str = r2_fastq}
  
  Array[Pair[File, File]] fastq_pairs = zip(fastq1.out, fastq2.out)
  
  scatter (fastq_pair in fastq_pairs) { call align_replicates {input: r1_fastq = fastq_pair.left, r2_fastq = fastq_pair.right, n_bp_trim_read1 = n_bp_trim_read1, n_bp_trim_read2 = n_bp_trim_read2, genome_index = genome_index, multicore = multicore, monitoring_script = monitoring_script,  memory = memory, disks = disks, cpu = cpu, preemptible = preemptible} }

  call merge_replicates {input: bams = align_replicates.bam, reports = align_replicates.output_report, samplename = samplename, chrom_sizes = chrom_sizes, genome_index = genome_index, multicore = multicore, monitoring_script = monitoring_script, memory = memory, disks = disks, cpu = cpu, preemptible = preemptible}
  
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
  File r1_fastq
  File r2_fastq
  File genome_index
  File monitoring_script
  
  String samplename = basename(r1_fastq)
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
    mkdir bismark_index
    tar zxf ${genome_index} -C bismark_index
    
		ln -s ${r1_fastq} r1.fastq.gz
		ln -s ${r2_fastq} r2.fastq.gz
		
		trim_galore --paired --clip_R1 ${n_bp_trim_read1} --clip_R2 ${n_bp_trim_read2} r1.fastq.gz r2.fastq.gz
        	
		bismark --genome bismark_index --multicore ${multicore} -1 r1_val_1.fq.gz -2 r2_val_2.fq.gz
    mv *bismark_bt2_pe.bam ${samplename}.bam
    mv *bismark_bt2_PE_report.txt ${samplename}_report.txt
    
    bismark2report --alignment_report ${samplename}_report.txt --output ${samplename}_bismark_report.html  
  
  }
  
  output {
    File bam = "${samplename}.bam"
    File output_report = "${samplename}_report.txt"
    File bismark_report_html = "${samplename}_bismark_report.html"
    File monitoring_log = "monitoring.log"
  }
  
  	runtime {
  	continueOnReturnCode: false
		docker: "aryeelab/bismark:latest"
		memory: memory
    disks: disks
    cpu: cpu
    preemptible: preemptible
 }
  
}

task merge_replicates {
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
		/src/Bismark-0.18.2/deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
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
  
  output {
            File monitoring_log = "monitoring.log" 
            File output_bam = "${samplename}.sorted.bam"
            File output_covgz = "${samplename}.bismark.cov.gz"
            File output_report = "${samplename}_report.txt"
            File mbias_report = "${samplename}.M-bias.txt"
            File output_bigwig = "${samplename}.bw"
            File output_bai = "${samplename}.sorted.bai"
  }
  
  runtime {
    continueOnReturnCode: false
    docker: "aryeelab/bismark:v1.0"
    memory: memory
    disks: disks
    cpu: cpu
    preemptible: preemptible
  }
}