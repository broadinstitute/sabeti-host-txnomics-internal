version 1.0

task inputValidation{
  input {
     Array[File] r1_fastq
     Array[File] r2_fastq
     Array[File] r3_fastq

     String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-umi_tagger:0.0.1"
  }

  command {
    set -e

    for cur_file in ~{ sep=' ' r1_fastq } ~{ sep=' ' r2_fastq } ~{ sep=' ' r3_fastq }; do
    	if [ -s $cur_file ]; then
	   echo File $cur_file is empty 
    	   exit 1
    	fi
    done
  }

  output {

  }  
  
}


task umiTagger {
  input {
      File r1_fastq
      File r2_fastq
      File r3_fastq

      String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-umi_tagger:0.0.1"
  }

  String r1_out_name = "r1.fastq"
  String r3_out_name = "r3.fastq"
  Int cpu = 1
  Int disk = ceil(size(r1_fastq, "GiB") * 4 + size(r2_fastq, "GiB") * 4 + size(r3_fastq, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e

    python3 /opt/tools/umi_tagger.py --r1 ~{r1_fastq} --r2 ~{r2_fastq} --r3 ~{r3_fastq} --r1-output ~{r1_out_name} --r3-output ~{r3_out_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File r1_out_fastq = r1_out_name
    File r3_out_fastq = r3_out_name
  }
}


task StarAlign {
  input {
      File r1_fastq
      File r2_fastq
      File reference
      String reference_prefix = ""
      String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-star:0.0.1"

      # runtime values
      Int machine_mem_mb = 50000
      Int cpu = 1
      # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
      Int disk = ceil(size(r1_fastq, "GiB") * 4 + size(r2_fastq, "GiB") * 4 + 100)
      # by default request non preemptible machine to make sure the slow star alignment step completes
      Int preemptible = 0
  }

  command {
      set -e

      # prepare reference
      mkdir genome_reference
      tar -xf "${reference}" -C genome_reference
      #rm "${reference}"
      
      find genome_reference -type d

      STAR \
          --runMode alignReads \
          --runThreadN ${cpu} \
          --genomeDir genome_reference/~{reference_prefix} \
          --readFilesIn ~{r1_fastq} ~{r2_fastq} \
          --outSAMtype BAM Unsorted \
          --outSAMattributes All \
          --outSAMunmapped None \
          --readFilesType Fastx \
          --runRNGseed 777

      # Extract some key metrics
      grep "Number of input reads" Log.final.out | cut -f 2 -d '|' | tr -d '\t\n' > outval_inputreads.txt
      grep "Uniquely mapped reads number" Log.final.out | cut -f 2 -d '|' | tr -d '\t\n' > outval_uniquelyaligned.txt
      grep "Number of reads mapped to multiple loci" Log.final.out | cut -f 2 -d '|' | tr -d '\t\n' > outval_multimaps.txt
      }

      runtime {
          docker: docker
          memory: "${machine_mem_mb} MiB"
          disks: "local-disk ${disk} HDD"
          cpu: cpu
          preemptible: preemptible
      }

  output {
      File out_bam = "Aligned.out.bam"
      File alignment_log = "Log.final.out"
      Int input_reads = read_int("outval_inputreads.txt")
      Int uniquely_aligned_reads = read_int("outval_uniquelyaligned.txt")
      Int multimap_reads = read_int("outval_multimaps.txt")
   }
}

task filterMultimaps {
    input {
        File input_bam
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 2 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-hisat2:0.0.1"

    String output_bam_filename = "output.bam"

    command {
        set -e
        samtools view -b -F 0xF00 ~{input_bam} > ~{output_bam_filename}
    }

    runtime {
        docker: docker
        memory: "4 GiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_bam = output_bam_filename
    }    
}

task estimate_library_complexity {
    input {
        File input_bam
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 2 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"

    String output_metrics_file = "est_lib_compl_metrics.txt"

    command {
        set -e
        java -jar "$PICARD_JAR_PATH" EstimateLibraryComplexity I=~{input_bam} O=~{output_metrics_file}
    }

    runtime {
        docker: docker
        memory: "8 GiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File library_complexity_metrics = output_metrics_file
    }
}

task countBamReads {
    input {
      File input_bam
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 1.1 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-hisat2:0.0.1"

    command {
      set -e
      samtools view -c ~{input_bam} > read_count.txt
    }

    runtime {
      docker: docker
      memory: "8 GiB"
      disks: "local-disk ${disk} HDD"
      cpu: cpu
      preemptible: preemptible
    }

    output {
      Int read_count = read_int("read_count.txt")
    }
}

task downsampleBam {
    input {
        File input_bam
        Float probability = 0.5
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 2 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"

    String downsampled_bam_file = "downsampled.bam"

    command {
        set -e
        java -Xms8G -Xmx14G -jar "$PICARD_JAR_PATH" DownsampleSam I=~{input_bam} O=~{downsampled_bam_file} STRATEGY=Chained P=~{probability}
    }

    runtime {
        docker: docker
        memory: "16 GiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_bam = downsampled_bam_file
    }
}

task sort_and_index {
    input {
       File input_bam
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-hisat2:0.0.1"

    String out_bam = "sorted.bam"

    command {
       set -e
       samtools sort ~{input_bam} > ~{out_bam}
       samtools index ~{out_bam}
    }

    runtime {
       docker: docker
       memory: "4 GiB"
       disks: "local-disk ~{disk} HDD"
       cpu: cpu
       preemptible: preemptible
    }

    output {
       File output_bam = out_bam
       File output_bam_index = out_bam + ".bai"
    }
}

task removeDuplicates {
    input {
       File bam
       File bam_index
    }

    Int cpu = 1
    Int disk = ceil(size(bam, "GiB") * 6 + 20)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-umi_tools:0.0.1"

    String out_bam = "deduplicated.bam"

    command {
       set -e
       #touch ~{bam}
       #sleep 1
       #touch ~{bam}.bai
       umi_tools dedup -I ~{bam} --chimeric-pairs=discard --unpaired-reads=discard --paired --output-stats=dedup_stats -S ~{out_bam}
    }

    runtime {
       docker: docker
       memory: "64 GiB"
       disks: "local-disk ~{disk} HDD"
       cpu: cpu
       preemptible: preemptible
    }

    output {
       File deduplicated_bam = out_bam
       Array[File] dedup_stats = glob("dedup_stats*")
       File dedup_stats_per_umi_per_position = "dedup_stats_per_umi_per_position.tsv"
    }

}

task generate_saturation_info {
    input {
        File dedup_stats_per_umi_per_position
	String sample_name = "Sample"
    }

    Int cpu = 1
    Int disk = 10
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-saturation:0.0.1"

    command {
        estimate_saturation.R -i ~{dedup_stats_per_umi_per_position} -n ~{sample_name}
    }

    runtime {
      docker: docker
      memory: "4 GiB"
      disks: "local-disk ~{disk} HDD"
      cpu: cpu
      preemptible: preemptible
    }

    output {
      File saturation_plot = "saturation_curve.pdf"
      File saturation_table = "saturation.tsv"
    }
    
}

task featureCounts {
    input {
        File input_bam
        File annotation_gtf
    }

    String counts_file_name = "counts.txt"

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 2 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/subread:2.0.1"

    command {
        featureCounts -p -T ~{cpu} -t exon -g gene_id -a ~{annotation_gtf} -o ~{counts_file_name} ~{input_bam}
    }

    runtime {
        docker: docker
        memory: "4 GiB"
        disks: "local-disk ~{disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File counts_file = counts_file_name
    }

}

task fastqMerge {
    input {
         Array[File] r1_fastq
         Array[File] r2_fastq
         Array[File] r3_fastq
    }
    
    Int cpu = 1
    Int disk = ceil(size(r1_fastq, "GiB") * 2.5 + size(r2_fastq, "GiB") * 2.5 + size(r3_fastq, "GiB") + 50)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-umi_tools:0.0.1"
    
    command {
       zcat ${sep=' ' r1_fastq} | gzip -c > r1_out.fastq.gz
       zcat ${sep=' ' r2_fastq} | gzip -c > r2_out.fastq.gz
       zcat ${sep=' ' r3_fastq} | gzip -c > r3_out.fastq.gz
    }
    
    runtime {
       docker: docker
       memory: "4 GiB"
       disks: "local-disk ~{disk} HDD"
       cpu: cpu
       preemptible: preemptible
    }

    output {
      File r1_fastq_out = "r1_out.fastq.gz"
      File r2_fastq_out = "r2_out.fastq.gz"
      File r3_fastq_out = "r3_out.fastq.gz"
    }

}

workflow umiRnaSeq {
  input {
      Array[File] r1_fastq
      Array[File] r2_fastq
      Array[File] r3_fastq
      File reference
      String reference_prefix
      File annotation_gtf
      
      Boolean subsampling = false
      Boolean do_estimate_library_complexity = false
      Boolean do_generate_saturation_info = true
  }

  Array[Float] subsample_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  String version = "umiRNASeq_v0.0.1"

  call inputValidation {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      r3_fastq = r3_fastq
  }

  call fastqMerge {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      r3_fastq = r3_fastq
  }

  call umiTagger {
    input:
      r1_fastq = fastqMerge.r1_fastq_out,
      r2_fastq = fastqMerge.r2_fastq_out,
      r3_fastq = fastqMerge.r3_fastq_out
  }

  call StarAlign {
    input:
      r1_fastq = umiTagger.r1_out_fastq,
      r2_fastq = umiTagger.r3_out_fastq,
      reference = reference,
      reference_prefix = reference_prefix
  }

  call filterMultimaps {
    input:
      input_bam = StarAlign.out_bam
  }

  if (do_estimate_library_complexity) {
    call estimate_library_complexity {
      input:
        input_bam = filterMultimaps.output_bam
    }
  }

  if ( subsampling ) {
    scatter ( p in subsample_values ) {
      call downsampleBam {
        input:
          input_bam = filterMultimaps.output_bam,
          probability = p
      }
      
      call countBamReads as preDedupCounting {
        input:
          input_bam = downsampleBam.output_bam
      }

      call sort_and_index as sortIndexSubsampled {
        input:
          input_bam = downsampleBam.output_bam
      }

      call removeDuplicates as removeDuplicatesSubsampled {
        input:
          bam = sortIndexSubsampled.output_bam,
          bam_index = sortIndexSubsampled.output_bam_index
      }

      call countBamReads {
        input:
          input_bam = removeDuplicatesSubsampled.deduplicated_bam
      }
    }
  }

  call sort_and_index {
    input:
      input_bam = filterMultimaps.output_bam
  }

  call removeDuplicates {
    input:
      bam = sort_and_index.output_bam,
      bam_index = sort_and_index.output_bam_index
  }

  if ( do_generate_saturation_info ) {
    call generate_saturation_info {
      input:
        dedup_stats_per_umi_per_position = removeDuplicates.dedup_stats_per_umi_per_position
    }
  }

  call featureCounts {
    input:
        input_bam = removeDuplicates.deduplicated_bam,
        annotation_gtf = annotation_gtf
  }

  output {
    String pipeline_version = version

    File r1_fastq_out = umiTagger.r1_out_fastq
    File r3_fastq_out = umiTagger.r3_out_fastq

    File aligned_bam = filterMultimaps.output_bam

    File? library_complexity_metrics = estimate_library_complexity.library_complexity_metrics

    File sorted_bam = sort_and_index.output_bam
    File sorted_bam_index = sort_and_index.output_bam_index

    File deduplicated_bam = removeDuplicates.deduplicated_bam
    Array[File] dedup_stats = removeDuplicates.dedup_stats
    File? saturation_plot = generate_saturation_info.saturation_plot
    File? saturation_table = generate_saturation_info.saturation_table

    Array[Float] subsample_values_out = subsample_values
    Array[Int]? pre_subsample_read_count = preDedupCounting.read_count
    Array[Int]? read_count = countBamReads.read_count
    
    File counts_file = featureCounts.counts_file

    Int alignment_input_reads = StarAlign.input_reads
    Int alignment_unique_reads = StarAlign.uniquely_aligned_reads
    Int alignment_multimaped_reads = StarAlign.multimap_reads
  }

}
