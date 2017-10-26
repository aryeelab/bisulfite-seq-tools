workflow call_bismark {
	call step1_bismark_wgbs { }
}

task step1_bismark_wgbs {
	File r1_fastq
	File r2_fastq
	String samplename
	File genome_index
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
            
		samtools sort -n -f ${samplename}.bam ${samplename}.sorted_by_readname.bam
		/src/Bismark-0.18.2/deduplicate_bismark -p --bam ${samplename}.sorted_by_readname.bam
		rm ${samplename}.sorted_by_readname.bam ${samplename}.bam
		mv ${samplename}.sorted_by_readname.deduplicated.bam ${samplename}.bam
            
		bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${samplename}.bam
		bismark2report --alignment_report ${samplename}_report.txt --output ${samplename}_bismark_report.html   
	}
        
	output {
		File output_covgz = "${samplename}.bismark.cov.gz"
		File output_report = "${samplename}_report.txt"
		File mbias_report = "${samplename}.M-bias.txt"
		File bismark_report_html = "${samplename}_bismark_report.html"
	}
        
	runtime {
		continueOnReturnCode: false
		docker: "sowmyaiyer/bismark_image:latest"
		memory: memory
		disks: disks
		cpu: cpu
		preemptible: preemptible
	}
}

