task combine_rda {
     Array[File] individual_rdas
     String set_name
     command {
             Rscript /Rscripts/combine_rda_wrapper.R --rdaFiles ${sep=',' individual_rdas} --outFile ${set_name}
     }
     runtime {
             docker: "aryeelab/methy"
     }
     output  {
             File combined_rda = "${set_name}.combined.rda"
     }
}
task qc_report {
     File combined_rda
     command {
             Rscript /Rscripts/generate_qc_report.R -f ${combined_rda}  -o $PWD 
     }
     runtime {
             docker: "aryeelab/methy"
     }
     output {
            File qc_report = "qcReport.html"
            }
}
workflow methpipeagg {
         String set_name
         call combine_rda {input: set_name = set_name}
         call qc_report {input: combined_rda = combine_rda.combined_rda}
}