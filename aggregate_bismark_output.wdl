workflow aggregate_bismark_output {
        Array[File] in_pe_reports_files
        Array[File] in_covgz_files
        Array[File] in_mbias_files
        File BSGenome_targz
        String BSGenome_package
        String Genome_build
                
        String memory
  		  String disks
  		  Int preemptible


        call step2_create_combined_bsseq {input:in_pe_reports_files=in_pe_reports_files,in_covgz_files=in_covgz_files,in_mbias_files=in_mbias_files,BSGenome_targz=BSGenome_targz,BSGenome_package=BSGenome_package,Genome_build=Genome_build, memory = memory, disks = disks, preemptible = preemptible}
}



task step2_create_combined_bsseq {
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
                tar -cf outputdir_final.tar outputdir_final
        }
        output {
                File final = "outputdir_final.tar"
                File final_qc_html = "qcReport.html"
        }
        runtime {
        continueOnReturnCode: false
        docker: "aryeelab/bsseq_aggregation:latest"
        memory: memory
        disks: disks
        preemptible: preemptible 
        }
}
