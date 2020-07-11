version 1.0

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
      Int disk = 200
      # by default request non preemptible machine to make sure the slow star alignment step completes
      Int preemptible = 0
  }

  command {
      set -e

      # prepare reference
      mkdir genome_reference
      tar -xf "${reference}" -C genome_reference
      rm "${reference}"
      
      find genome_reference -type d

      STAR \
          --runMode alignReads \
          --runThreadN ${cpu} \
          --genomeDir genome_reference/~{reference_prefix} \
          --readFilesIn ~{r1_fastq} ~{r2_fastq} \
          --outSAMtype BAM Unsorted \
          --outSAMattributes All \
          --outSAMunmapped Within \
          --readFilesType Fastx \
          --runRNGseed 777
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
    Int disk = ceil(size(bam, "GiB") * 4 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-umi_tools:0.0.1"

    String out_bam = "deduplicated.bam"

    command {
       set -e
       touch ~{bam}
       sleep 1
       touch ~{bam}.bai
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
    }

}

task fastqMerge {
    input {
         Array[File] r1_fastq
         Array[File] r2_fastq
         Array[File] r3_fastq
    }
    
    Int cpu = 1
    Int disk = 500
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
      File r1_fastq = "r1_out.fastq.gz"
      File r2_fastq = "r2_out.fastq.gz"
      File r3_fastq = "r3_out.fastq.gz"
    }

}

workflow umiRnaSeq {
  input {
      Array[File] r1_fastq
      Array[File] r2_fastq
      Array[File] r3_fastq
      File reference
  }

  String version = "umiRNASeq_v0.0.1"
  
  call fastqMerge {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      r3_fastq = r3_fastq
  }

  call umiTagger {
    input:
      r1_fastq = fastqMerge.r1_fastq,
      r2_fastq = fastqMerge.r2_fastq,
      r3_fastq = fastqMerge.r3_fastq
  }

  call StarAlign {
    input:
      r1_fastq = umiTagger.r1_out_fastq,
      r2_fastq = umiTagger.r3_out_fastq,
      reference = reference
  }

  call sort_and_index {
    input:
      input_bam = StarAlign.out_bam
  }

  call removeDuplicates {
    input:
      bam = sort_and_index.output_bam,
      bam_index = sort_and_index.output_bam_index
  }

  output {
    String pipeline_version = version

    File r1_fastq_out = umiTagger.r1_out_fastq
    File r3_fastq_out = umiTagger.r3_out_fastq

    File aligned_bam = StarAlign.out_bam

    File sorted_bam = sort_and_index.output_bam
    File sorted_bam_index = sort_and_index.output_bam_index

    File deduplicated_bam = removeDuplicates.deduplicated_bam
    Array[File] dedup_stats = removeDuplicates.dedup_stats
  }

}
