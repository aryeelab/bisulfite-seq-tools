task bismark {
  File fastq1
  File fastq2
  File genome_index
  String sample_id
  File monitoring_script

  command {    
  		chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
        mkdir bismark_index
        tar zxf ${genome_index} -C bismark_index
        bismark --genome bismark_index --basename ${sample_id} -1 ${fastq1} -2 ${fastq2}
        bismark_methylation_extractor --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${sample_id}_pe.bam
  }
  runtime {
          docker: "aryeelab/bismark"
          memory: "24 GB"
          disks: "local-disk 100 SSD"
          cpu: 4
          preemptible: 1
  }
  output {
         File cov = "${sample_id}_pe.bismark.cov.gz"
         File report = "${sample_id}_PE_report.txt"
         File monitoring_log = "monitoring.log" 
  }

}


workflow bspipe {
    File fastq1
    File fastq2
    File genome_index
    String sample_id
    File monitoring_script
    call bismark {input: fastq1 = fastq1, fastq2 = fastq2, sample_id = sample_id, genome_index = genome_index, monitoring_script = monitoring_script}
}
