workflow aggregate_bismark_output {
	Array[File] in_pe_reports_files
	Array[File] in_covgz_files
	String BSGenome_targz
	String BSGenome_package
	
	call step2_create_combined_bsseq {input:in_pe_reports_files=in_pe_reports_files,in_covgz_files=in_covgz_files,genome=genome,BSGenome_targz=BSGenome_targz,BSGenome_package=BSGenome_package}
}



task step2_create_combined_bsseq {
	Array[File] in_pe_reports_files
	Array[File] in_covgz_files
	String BSGenome_targz
	String BSGenome_package
	command {
		R CMD INSTALL BSGenome_targz
		mkdir outputdir_final
		echo """${sep="\n" in_pe_reports_files}""" > all_pe_reports.txt
		echo """${sep="\n" in_covgz_files}""" > all_cov_gz.txt
		Rscript --vanilla /bismark_to_bsseq_file_list.R -i all_pe_reports.txt -j all_cov_gz.txt --bsgenome ${BSGenome_package} -o outputdir_final
		Rscript --vanilla /qc_report.R outputdir_final ${BSGenome_package}
		tar -cf outputdir_final.tar outputdir_final
	}
	output {
		File final = "outputdir_final.tar"
		File final_qc_html = "qcReport.html"
	}
	runtime {
	continueOnReturnCode: false
	docker: "sowmyaiyer/methyl_aggregation:latest"
	memory: "16GB"
	}
}

task qc_reports {
	File directory
	String genome
	command {
	mkdir untarred_directory
	tar xvf ${directory} -C untarred_directory
	Rscript -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"
	Rscript --vanilla /qc_report.R untarred_directory ${genome}
	}
	output {
		File final = "qcReport.html"
	}
	runtime {
	continueOnReturnCode: false
	docker: "sowmyaiyer/methyl_aggregation:latest"
	memory: "16GB"
	}
}
