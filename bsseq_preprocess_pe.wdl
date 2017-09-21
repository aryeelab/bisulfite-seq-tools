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
        bismark --genome bismark_index -1 ${fastq1} -2 ${fastq2}
        # The file renaming is necessary since this version of bismark doesn't allow the 
        # use of --multicore with --basename
        mv *bismark_bt2_pe.bam ${sample_id}.bam
        mv *bismark_bt2_PE_report.txt ${sample_id}_report.txt
        bismark_methylation_extractor --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${sample_id}.bam
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
         File report = "${sample_id}_report.txt"
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
