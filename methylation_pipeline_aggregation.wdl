task combine_read_metrics{
     Array[File] individual_read_metrics
     command{
	     cat ${sep=' ' individual_read_metrics} > read_metrics.txt
     }
     runtime {
	     docker: "adunford/methy:3"
     }
     output  {
     	     File combined_metrics = "read_metrics.txt"
     }
}

task combine_rda {
     Array[File] individual_rdas
     String set_name

     command {
     	     Rscript /RScripts/combine_rda.R --rdaFiles ${sep=',' individual_rdas} --outFile ${set_name}
     }
     runtime {
             docker: "adunford/methy:3"
     }
     output  {
     	     File combined_rda = "${set_name}.combined.rda"
     }
}


workflow methpipeagg {
	 String set_name
	 call combine_read_metrics
	 call combine_rda {input: set_name = set_name}
}