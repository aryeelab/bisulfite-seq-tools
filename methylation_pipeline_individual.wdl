task bsmap {
  File fastq1
  File fastq2
  File ref_genome
  String sample


  command {
         bsmap -a ${fastq1} -b ${fastq2} -d ${ref_genome} -p 3 -v 0.05 -s 16 -r 0 -u -S 1 -R -o ${sample}_raw_bs.bam
  }
  runtime {
          docker: "aryeelab/methy"
          memory: "16 GB"
          disks: "local-disk 100 SSD"
  }
  output {
         File raw_bs_bam = "${sample}_raw_bs.bam"
         #File genome = "${ref_genome}"
         #String sample_id = "${sample}"
  }
}

task samtools_sort {
  File raw_bs_bam
  String sample_id
  command {
          samtools sort ${raw_bs_bam} ${sample_id}_bs.sorted && samtools index ${sample_id}_bs.sorted.bam
  }
  runtime {
          docker: "aryeelab/methy"
          memory: "16 GB"
          disks: "local-disk 100 SSD"
  }
  output {
         File sorted_bs_bam   = "${sample_id}_bs.sorted.bam"
         File sorted_bs_bai   = "${sample_id}_bs.sorted.bam.bai"
	}
}

task samtools_read_metrics{
     File sorted_bs_bam
     String sample_id
     command{
        echo ${sample_id} `samtools view ${sorted_bs_bam} | wc -l` `samtools view -F 4 ${sorted_bs_bam} | wc -l` > ${sample_id}.read_metrics.txt
     }
     runtime{
        docker: "aryeelab/methy"
        memory: "16 GB"
        disks: "local-disk 100 SSD"
     }
     output {
        File read_metrics = "${sample_id}.read_metrics.txt"
     }
}

task MethylDackel {
        File genome
        File sorted_bs_bam
        String sample_id
        command {
                MethylDackel extract ${genome} ${sorted_bs_bam} -o ${sample_id}
                grep -v '^track' ${sample_id}_CpG.bedGraph  > tmp
                mv tmp ${sample_id}_CpG.bedGraph
        }
        runtime {
              docker: "aryeelab/methy"
    	      memory: "16 GB"
	          defaultDisks: "local-disk 100 SSD"

        }
        output {
                File bed = "${sample_id}_CpG.bedGraph"
        }
}

task MethylDackel_CHH {
     File genome
     File sorted_bs_bam
     String sample_id
     command {
             MethylDackel extract --CHH ${genome} ${sorted_bs_bam} -o ${sample_id}
     }
     runtime {
             docker: "aryeelab/methy"
             memory: "16 GB"
             disks: "local-disk 100 SSD"

     }
     output {
            File chh_bed = "${sample_id}_CHH.bedGraph"
     }
}
  String sample
task bs_conversion_rate{
     File chh_bed
     String sample_id
     command {
             sh /executable_files/collect_bsconv_metrics.sh ${sample_id} ${chh_bed} 
     }
     runtime {
             docker: "aryeelab/methy"
     }
     output{
             File bsconv = "${sample_id}_bsconv.txt"
     }
}
task create_rda {
     File bed
     String sample_id
     File bsconv
     File read_metrics
     command{
             Rscript /Rscripts/create_rda_wrapper.R -f ${bed} -o ${sample_id}.rda -b ${bsconv} -r ${read_metrics}
     }
     runtime {
             docker: "aryeelab/methy"
     }
     output {
            File rda = "${sample_id}.rda"
     }
}
workflow methpipeindv {
         File ref_genome
         File sample_id
         call bsmap     {input: sample = sample_id, ref_genome = ref_genome}
         call samtools_sort {input: raw_bs_bam = bsmap.raw_bs_bam, sample_id = sample_id }
         call samtools_read_metrics {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample_id}
         call MethylDackel {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample_id, genome = ref_genome}
         call MethylDackel_CHH  {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample_id, genome = ref_genome}
         call bs_conversion_rate {input: chh_bed = MethylDackel_CHH.chh_bed, sample_id = sample_id}
         call create_rda    {input: bed = MethylDackel.bed, sample_id = sample_id,bsconv = bs_conversion_rate.bsconv, read_metrics = samtools_read_metrics.read_metrics}
}
