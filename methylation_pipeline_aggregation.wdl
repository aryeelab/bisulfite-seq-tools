task bsmap {
  File fastq1
  File fastq2
  File genome
  String sample


  command {
  	 bsmap -a ${fastq1} -b ${fastq2} -d ${genome} -p 4 -v 0.05 -s 16 -r 0 -u -S 1 -R -o ${sample}_raw_bs.bam
  }
  runtime {
  	  docker: "adunford/methy:9"
	  #memory: "16 GB"
	  #defaultDisks: "local-disk 100 SSD"
  }
  output {
	 File raw_bs_bam = "./${sample}_raw_bs.bam"


  }
}

task samtools_sort {
  File raw_bs_bam
  String sample_id
  command {
          samtools sort ${raw_bs_bam} ${sample_id}_bs.sorted && samtools index ${sample_id}_bs.sorted.bam
  }
  runtime {
          docker: "adunford/methy:9"
          #memory: "16 GB"
          #defaultDisks: "local-disk 100 SSD"
  }
  output {
         File sorted_bs_bam   = "${sample_id}_bs.sorted.bam"

  }
}

task samtools_read_metrics{
     File sorted_bs_bam
     String sample_id
     command{
	echo ${sample_id} `samtools view ${sorted_bs_bam} | wc -l` `samtools view -F 4 ${sorted_bs_bam} | wc -l` > ${sample_id}.read_metrics.txt
     }
     runtime{
	docker: "adunford/methy:9"
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
		docker: "adunford/methy:9"
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
     	     docker: "adunford/methy:9"
     }
     output {
     	    File chh_bed = "${sample_id}_CHH.bedGraph"
     }
}

task bs_conversion_rate{
     File chh_bed
     String sample_id
     command {
	     sh /executable_files/collect_bsconv_metrics.sh ${sample_id} ${chh_bed}
     }
     runtime {
     	     docker: "adunford/methy:9"
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
     	     docker: "adunford/methy:9"
     }
     output {
     	    File rda = "${sample_id}.rda"
     }

}

task combine_rda {
     Array[File] individual_rdas
     String set_name

     command {
                Rscript /Rscripts/combine_rda_wrapper.R --rdaFiles ${sep=',' individual_rdas} --outFile ${set_name}
	}
	runtime {
		docker: "adunford/methy:9"
	}
	output  {
		File combined_rda = "${set_name}.combined.rda"
	}
}

task qc_report {
     File combined_rda
     #String organism
     #String outdir
     #File genome
     command {
     	     Rscript /Rscripts/generate_qc_report.R -f ${combined_rda}  -o $PWD 
     }
     runtime {
	     docker: "adunford/methy:9"
     }
     output {
     	    File qc_report = "qcReport.html"
     }
}

workflow methpipeagg {
	 File genome
	 String set_name
	 Array[Map[String,String]] samples
	 scatter (sample in samples){
	 	 call bsmap	{input: sample = sample["sample_id"], fastq1 = sample["fastq1"], fastq2 = sample["fastq2"], genome = genome}
	 	 call samtools_sort {input: raw_bs_bam = bsmap.raw_bs_bam, sample_id = sample["sample_id"] }
	 	 call samtools_read_metrics {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample["sample_id"]}
	 	 call MethylDackel {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample["sample_id"], genome = genome}
	 	 call MethylDackel_CHH  {input: sorted_bs_bam = samtools_sort.sorted_bs_bam, sample_id = sample["sample_id"], genome = genome}
	  	 call bs_conversion_rate {input: chh_bed = MethylDackel_CHH.chh_bed, sample_id = sample["sample_id"]}
	  	 call create_rda    {input: bed = MethylDackel.bed, sample_id = sample["sample_id"], bsconv = bs_conversion_rate.bsconv, read_metrics = samtools_read_metrics.read_metrics}
	 }
	 call combine_rda {input: individual_rdas = create_rda.rda, set_name = set_name}
	 call qc_report {input: combined_rda =  combine_rda.combined_rda}
	 
}