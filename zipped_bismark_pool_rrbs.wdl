workflow call_bismark_pool {

  String r1_fastq
  String r2_fastq
  File genome_index
  File monitoring_script
  
  String multicore
  
  
  
   # Split the comma-separated string of fastq file names into arrays
  call split_string_into_array as fastq1 {input: str = r1_fastq}
  call split_string_into_array as fastq2 {input: str = r2_fastq}
  
  Array[Pair[File, File]] fastq_pairs = zip(fastq1.out, fastq2.out)
  
  scatter (fastq_pair in fastq_pairs) { call align_replicates {input: r1_fastq = fastq_pair.left, r2_fastq = fastq_pair.right, genome_index = genome_index, multicore = multicore, monitoring_script = monitoring_script} }
  
}

task split_string_into_array {
    String str 
    String arr = "{ARR[@]}"
    command {
        IFS=',' read -ra ARR <<< "${str}"
        for i in "$${arr}"; do
            echo "$i" | tr -d " " >> out
        done
    }

    runtime {
        docker: "debian:stretch"
    }
    
    output {
        Array[String] out = read_lines("out")
    }
}

task align_replicates{
  File r1_fastq
  File r2_fastq
  File genome_index
  File monitoring_script
  
  String multicore

  
  command {
    chmod u+x ${monitoring_script}
    ${monitoring_script} > monitoring.log &
    mkdir bismark_index
    tar zxf ${genome_index} -C bismark_index
    
    bismark --genome bismark_index --multicore ${multicore} -1 ${r1_fastq} -2 ${r2_fastq}
  
  }
  
  output {
    File bam = "*bismark_bt2_pe.bam.bam"
    File monitoring_log = "monitoring.log"
  }
  
}


