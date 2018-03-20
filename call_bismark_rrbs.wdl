workflow call_bismark {
        call step1_bismark_rrbs { }
    }

task step1_bismark_rrbs {
        File r1_fastq
        File r2_fastq
        String samplename
        File genome_index
        File chrom_sizes
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
            
            ln -s ${r1_fastq} r1.fastq.gz
		        ln -s ${r2_fastq} r2.fastq.gz
            
            trim_galore --rrbs --paired r1.fastq.gz r2.fastq.gz
            
            bismark --genome bismark_index --multicore ${multicore} -1 r1_val_1.fq.gz -2 r2_val_2.fq.gz
            # The file renaming below is necessary since this version of bismark doesn't allow the 
            # use of --multicore with --basename

            
            mv *bismark_bt2_pe.bam ${samplename}.bam
            mv *bismark_bt2_PE_report.txt ${samplename}_report.txt
            
            
  
            samtools sort -o ${samplename}.sorted.bam ${samplename}.bam
            samtools index ${samplename}.sorted.bam ${samplename}.sorted.bai
          
            
            
            bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${samplename}.bam
            bismark2report --alignment_report ${samplename}_report.txt --output ${samplename}_bismark_report.html  
            gunzip "${samplename}.bedGraph.gz"
            
            
            #LC_COLLATE=C sort -k1,1 -k2,2n ${samplename}.bedGraph > ${samplename}.sorted.bedGraph
            bedtools sort -i ${samplename}.bedGraph > ${samplename}.sorted.bedGraph
            /usr/local/bin/bedGraphToBigWig ${samplename}.sorted.bedGraph ${chrom_sizes} ${samplename}.bw
        }
                
        output {
            File output_bam = "${samplename}.sorted.bam"
            File output_covgz = "${samplename}.bismark.cov.gz"
            File output_report = "${samplename}_report.txt"
            File mbias_report = "${samplename}.M-bias.txt"
            File output_bigwig = "${samplename}.bw"
            File bismark_report_html = "${samplename}_bismark_report.html"
            File output_bai = "${samplename}.sorted.bai"
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
