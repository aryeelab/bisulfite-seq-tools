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
    samtools sort -n -o ${samplename}.sorted_by_readname.bam ${samplename}.merged.bam
		/src/Bismark-0.18.2/deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
		rm ${samplename}.sorted_by_readname.bam ${samplename}.bam
		mv ${samplename}.sorted_by_readname.deduplicated.bam ${samplename}.sorted.bam
    
    
  }
  
  output {
    File output_bam = "${samplename}.sorted.bam"
  }
  
  runtime {
    continueOnReturnCode: false
    docker: "aryeelab/bismark:latest"
    memory: "20GB"
    disks: "local-disk 500 SSD" 
  }
}