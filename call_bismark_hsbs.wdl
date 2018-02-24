workflow call_bismark {
        	        call step1_bismark_hsbs { }
}
task step1_bismark_hsbs {
        File r1_fastq
        File r2_fastq
        String samplename
        File genome_index
        File target_region_bed
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
            bismark --genome bismark_index --multicore ${multicore} -1 ${r1_fastq} -2 ${r2_fastq}
            # The file renaming below is necessary since this version of bismark doesn't allow the 
            # use of --multicore with --basename
            mv *bismark_bt2_pe.bam ${samplename}.bam
            mv *bismark_bt2_PE_report.txt ${samplename}_report.txt
            
            samtools sort -n -o ${samplename}.sorted_by_readname.bam ${samplename}.bam 
            /src/Bismark-0.18.2/deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
            samtools sort -o ${samplename}.sorted.deduplicated.bam ${samplename}.sorted_by_readname.deduplicated.bam 
            mv ${samplename}.sorted.deduplicated.bam ${samplename}.bam
            samtools index ${samplename}.bam ${samplename}.bai 
            
            bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${samplename}.bam
            bismark2report --alignment_report ${samplename}_report.txt --output ${samplename}_bismark_report.html

            # Target coverage report            
            # Total reads?
            TOTAL_READS=$(samtools view -F 4 ${samplename}.bam | wc -l | tr -d '[:space:]')
            echo "# total_reads=$TOTAL_READS" | tee ${samplename}_target_coverage.bed
            # How many reads overlap targets?
            READS_OVERLAPPING_TARGETS=$(bedtools intersect -u -bed -a ${samplename}.bam -b ${target_region_bed} | wc -l | tr -d '[:space:]')
            echo "# reads_overlapping_targets=$READS_OVERLAPPING_TARGETS" | tee -a ${samplename}_target_coverage.bed
            # How many reads per target?
            bedtools intersect -c -a ${target_region_bed} -b ${samplename}.bam >> ${samplename}_target_coverage.bed

        }
        output {
            File output_covgz = "${samplename}.bismark.cov.gz"
            File output_report = "${samplename}_report.txt"
            File output_bam = "${samplename}.bam"
            File output_bai = "${samplename}.bai"
		        File mbias_report = "${samplename}.M-bias.txt"
		        File bismark_report_html = "${samplename}_bismark_report.html"
            File target_coverage_report = "${samplename}_target_coverage.bed"
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

