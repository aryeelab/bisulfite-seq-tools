task bismark {
  File fastq
  File genome_index
  String sample_id
  File monitoring_script
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
        bismark --genome bismark_index --multicore ${multicore} ${fastq}
        # The file renaming is necessary since this version of bismark doesn't allow the 
        # use of --multicore with --basename
        mv *_bismark_bt2.bam ${sample_id}.bam
        mv *bismark_bt2_SE_report.txt ${sample_id}_report.txt
        bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${sample_id}.bam
  }
  runtime {
          docker: "aryeelab/bismark"
          memory: memory
          disks: disks
          cpu: cpu
          preemptible: preemptible
  }
  output {
         File cov = "${sample_id}.bismark.cov.gz"
         File report = "${sample_id}_report.txt"
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
