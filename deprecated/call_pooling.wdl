workflow call_bismark_pool {
  Array[File] in_bam_files
  String samplename

        call merge_replicates {input:in_bam_files=in_bam_files,samplename=samplename}
}

task merge_replicates {
  Array[File] in_bam_files
  String samplename

  command {
    
    samtools merge ${samplename}.merged.bam ${sep=' ' in_bam_files}
    
  
  }
  
  output {
    File output_bam = "${samplename}.merged.bam"
  }
  
  runtime {
        continueOnReturnCode: false
        docker: "aryeelab/bismark:latest"
        memory: "20GB"
        disks: "local-disk 500 SSD" 
  }
}