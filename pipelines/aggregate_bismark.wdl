workflow aggregate_bismark_output {

    String version = "dev"

    # 'dev' pipeline versions use the image with the 'latest' tag.
    # release pipeline versions use images tagged with the version number itself
    String image_id = sub(version, "dev", "latest")

    Array[File] in_pe_reports_files
    Array[File] in_covgz_files
    Array[File] in_mbias_files
    File BSGenome_targz
    String BSGenome_package
    String Genome_build

    String memory
    String disks
    Int preemptible

    call version_info { input: image_id = image_id }

    call create_combined_bsseq {input:  in_pe_reports_files=in_pe_reports_files,
                                    in_covgz_files=in_covgz_files,
                                    in_mbias_files=in_mbias_files,
                                    BSGenome_targz=BSGenome_targz,
                                    BSGenome_package=BSGenome_package,
                                    Genome_build=Genome_build,
                                    memory = memory, disks = disks, preemptible = preemptible, image_id = image_id}
    }

task version_info {
	String image_id
	command <<<
		cat /VERSION
	>>>
	runtime {
            continueOnReturnCode: false
            docker: "gcr.io/aryeelab/bismark:${image_id}"
            cpu: 1
            memory: "1GB"
        }
	output {
	    String pipeline_version = read_string(stdout())
        }
}


task create_combined_bsseq {
        String image_id
        Array[File] in_pe_reports_files
        Array[File] in_covgz_files
        Array[File] in_mbias_files
        File BSGenome_targz
        String BSGenome_package
        String Genome_build
        
        String memory
  		String disks
  		Int preemptible
  		  
        command {
                R CMD INSTALL ${BSGenome_targz}
                mkdir outputdir_final
                mkdir mbias_files
                echo """${sep="\n" in_pe_reports_files}""" > all_pe_reports.txt
                echo """${sep="\n" in_covgz_files}""" > all_cov_gz.txt
                echo """${sep="\n" in_mbias_files}""" > all_mbias.txt
                Rscript --vanilla /bismark_to_bsseq_file_list.R -i all_pe_reports.txt -j all_cov_gz.txt -k all_mbias.txt --bsgenome ${BSGenome_package} -o outputdir_final -m mbias_files
                Rscript -e "library(dplyr);library(scmeth);library('${BSGenome_package}', character.only = TRUE);bs <- HDF5Array::loadHDF5SummarizedExperiment(dir='outputdir_final');getwd();report(bs, '/', get('${BSGenome_package}'), '${Genome_build}',mbiasDir='mbias_files')"
                mv /qcReport.html .
                mv /QC_Summary.txt .
                mv /Downsample_Table.txt .
                mv /mbias_Table.txt .
                tar -cf outputdir_final.tar outputdir_final
        }
        output {
                File final = "outputdir_final.tar"
                File final_qc_html = "qcReport.html"
                File final_qc_table = "QC_Summary.txt"
                File final_downsample_table = "Downsample_Table.txt"
                File final_mbias_table = "mbias_Table.txt"
        }
        runtime {
        continueOnReturnCode: false
    	docker: "gcr.io/aryeelab/bismark:${image_id}"
        memory: memory
        disks: disks
        preemptible: preemptible 
        }
}
