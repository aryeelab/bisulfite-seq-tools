workflow aggregate_bismark_output {
        
        call step2_create_combined_bsseq {  }
}



task step2_create_combined_bsseq {
	Array[File] in_pe_reports_files
	Array[File] in_covgz_files
	File genome_CG_Robj
	String BSgenome_name
        command {
		mkdir outputdir_final
		echo """${sep="\n" in_pe_reports_files}""" > all_pe_reports.txt
                echo """${sep="\n" in_covgz_files}""" > all_cov_gz.txt
		Rscript --vanilla /bismark_to_bsseq_file_list.R -i all_pe_reports.txt -j all_cov_gz.txt --genome ${genome_CG_Robj} -o outputdir_final
		tar -cf outputdir_final.tar outputdir_final
        }
        output {
                File final = "outputdir_final.tar"
        }
        runtime {
	continueOnReturnCode: false
        docker: "sowmyaiyer/methyl_aggregation"
        memory: "16GB"
        }
}
