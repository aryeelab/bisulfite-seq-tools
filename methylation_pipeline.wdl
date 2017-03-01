task bsmap {
  File fastq1
  File fastq2
  File ref_genome
  String sample


  command {
  	 bsmap -a ${fastq1} -b ${fastq2} -d ${ref_genome} -p 4 -v 0.05 -s 16 -r 0 -u -S 1 -R -o ${sample}_raw_bs.bam
  }
  runtime {
  	  docker: "adunford/methy:1"
	  #memory: "16 GB"
	  #defaultDisks: "local-disk 100 SSD"
  }
  output {
	 File raw_bs_bam = "./${sample}_raw_bs.bam"
	 File genome = "${ref_genome}"
	 String sample_id = "${sample}"
  }
}

task samtools_sort {
  File raw_bs_bam
  String sample_id
  command {
          samtools sort ${raw_bs_bam} ${sample_id}_bs.sorted && samtools index ${sample_id}_bs.sorted.bam
  }
  runtime {
          docker: "adunford/methy:1"
          #memory: "16 GB"
          #defaultDisks: "local-disk 100 SSD"
  }
  output {
         File sorted_bs_bam   = "${sample_id}_bs.sorted.bam"
	 #File sorted_bs_index = "${sample_id}_bs.sorted.bam.bai"
  }
}


task MethylDackel {
     File genome
     File sorted_bs_bam
     String sample_id
     command {
     	     MethylDackel extract ${genome} ${sorted_bs_bam} -o ${sample_id}
     }
     runtime {
     	     docker: "adunford/methy:1"
     }
     output {
     	    File bed = "${sample_id}_CpG.bedGraph"
     }
}


workflow methpipe {
	 call bsmap
	 call samtools_sort {input: raw_bs_bam = bsmap.raw_bs_bam, sample_id = bsmap.sample_id }
	 call MethylDackel {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = bsmap.sample_id, genome = bsmap.genome}
}