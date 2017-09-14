task bismark {
  File fastq
  File genome_index
  String sample_id
  File monitoring_script

  command {    
  		chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
        mkdir bismark_index
        tar zxf ${genome_index} -C bismark_index
        bismark --genome bismark_index --basename ${sample_id} ${fastq} 
        bismark_methylation_extractor --gzip --bedGraph --buffer_size 4G --genome_folder bismark_index ${sample_id}.bam
  }
  runtime {
          docker: "aryeelab/bismark"
          memory: "24 GB"
          disks: "local-disk 100 SSD"
          cpu: 4
          preemptible: 1
  }
  output {
         File cov = "${sample_id}.bismark.cov.gz"
         File report = "${sample_id}_SE_report.txt"
         File monitoring_log = "monitoring.log" 
  }

}


workflow bspipe {
    File fastq
    File genome_index
    String sample_id
    File monitoring_script
    call bismark {input: fastq = fastq, sample_id = sample_id, genome_index = genome_index, monitoring_script = monitoring_script}
}
