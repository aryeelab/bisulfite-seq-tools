workflow call_bismark {
        call step1_bismark_rrbs { }
    }

task step1_bismark_rrbs {
        File r1_fastq
        File r2_fastq
        String samplename
        File genome_index
        File monitoring_script
        
        command {
            chmod u+x ${monitoring_script}
            ${monitoring_script} > monitoring.log &
            mkdir bismark_index
            tar zxf ${genome_index} -C bismark_index
            bismark --genome bismark_index --basename  ${samplename} -1 ${r1_fastq} -2 ${r2_fastq}
            bismark_methylation_extractor --gzip --bedGraph --buffer_size 4G --genome_folder bismark_index ${samplename}_pe.bam
            bismark2report --output ${samplename}_bismark_report.html
        }
                
        output {
            File output_covgz = "${samplename}_pe.bismark.cov.gz"
            File output_pe_report = "${samplename}_PE_report.txt"
            File mbias_report = "${samplename}_pe.M-bias.txt"
            File bismark_report_html = "${samplename}_bismark_report.html"
        }
        
        runtime {
            continueOnReturnCode: false
            docker: "sowmyaiyer/methyl:latest"
        }
}
