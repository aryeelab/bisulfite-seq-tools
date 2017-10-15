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
        command {
		    chmod u+x ${monitoring_script}
	        ${monitoring_script} > monitoring.log &
        	mkdir bismark_index
        	tar zxf ${genome_index} -C bismark_index
            bismark --genome bismark_index --basename  ${samplename} -1 ${r1_fastq} -2 ${r2_fastq}
            samtools sort -n -f ${samplename}_pe.bam ${samplename}_pe.sorted_by_readname.bam
            /src/Bismark-0.18.2/deduplicate_bismark -p --bam ${samplename}_pe.sorted_by_readname.bam
            samtools sort -f ${samplename}_pe.sorted_by_readname.deduplicated.bam ${samplename}_pe.sorted.deduplicated.bam
            mv ${samplename}_pe.sorted.deduplicated.bam ${samplename}_pe.bam
            bismark_methylation_extractor --gzip --bedGraph --buffer_size 4G --genome_folder bismark_index ${samplename}_pe.bam
            bismark2report --output ${samplename}_bismark_report.html

            # Target coverage report            
            # Total reads?
            TOTAL_READS=$(samtools view -F 4 ${samplename}_pe.bam | wc -l | tr -d '[:space:]')
            echo "# total_reads=$TOTAL_READS" | tee ${samplename}_target_coverage.bed
            # How many reads overlap targets?
            READS_OVERLAPPING_TARGETS=$(bedtools intersect -u -bed -a ${samplename}_pe.bam -b ${target_region_bed} | wc -l | tr -d '[:space:]')
            echo "# reads_overlapping_targets=$READS_OVERLAPPING_TARGETS" | tee -a ${samplename}_target_coverage.bed
            # How many reads per target?
            bedtools intersect -c -a ${target_region_bed} -b ${samplename}_pe.bam >> ${samplename}_target_coverage.bed

        }
        output {
            File output_covgz = "${samplename}_pe.bismark.cov.gz"
            File output_pe_report = "${samplename}_PE_report.txt"
		    File mbias_report = "${samplename}_pe.M-bias.txt"
		    File bismark_report_html = "${samplename}_bismark_report.html"
            File target_coverage_report = "${samplename}_target_coverage.bed"
        }
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/bismark_image:latest"
        }
}

