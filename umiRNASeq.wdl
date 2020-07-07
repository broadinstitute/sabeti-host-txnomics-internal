version 1.0

task umiTagger {
  input {
      File r1_fastq
      File r2_fastq
      File r3_fastq

      String docker = "umi_tagger:0.0.1"
  }

  String r1_out_name = "r1.fastq.gz"
  String r3_out_name = "r3.fastq.gz"
  Int cpu = 1
  Int disk = ceil(size(r1_fastq, "GiB") * 2 + size(r2_fastq, "GiB") + size(r3_fastq, "GiB") * 2 +  1)
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

task hisat2_align_pe {
   input {
     File r1_fastq
     File r2_fastq
     File reference
     String reference_prefix = "grch38/genome"
     String docker = "hisat2:0.0.1"
   }

  Int cpu = 1
  Int disk = ceil(size(r1_fastq, "GiB") * 2 + size(r2_fastq, "GiB") + size(reference, "GiB") * 2 +  10)
  Int preemptible = 3


   String out_bam = "aligned.bam"

   command {
     set -e
     tar xzf ~{reference}
     
     hisat2 -x ~{reference_prefix} -1 ~{r1_fastq} -2 ~{r2_fastq} | samtools view -b > ~{out_bam}
   }

   runtime {
     docker: docker
     memory: "4 GiB"
     disks: "local-disk ~{disk} HDD"
     cpu: cpu
     preemptible: preemptible
   }

   output {
     File out_bam = out_bam
   }
}


task sort_and_index {
    input {
       File input_bam
    }

    Int cpu = 1
    Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
    Int preemptible = 3
    String docker = "hisat2:0.0.1"

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
    String docker = "umi_tools:0.0.1"

    String out_bam = "deduplicated.bam"

    command {
       set -e
       umi_tools dedup -I ~{bam} --paired --output-stats=dedup_stats -S ~{out_bam}
    }

    runtime {
       docker: docker
       memory: "4 GiB"
       disks: "local-disk ~{disk} HDD"
       cpu: cpu
       preemptible: preemptible
    }

    output {
       File deduplicated_bam = out_bam
       Array[File] dedup_stats = glob("dedup_stats*")
    }

}

workflow umiRnaSeq {
  input {
      File r1_fastq
      File r2_fastq
      File r3_fastq
      File reference
  }

  String version = "umiRNASeq_v0.0.1"

  call umiTagger {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      r3_fastq = r3_fastq
  }

  call hisat2_align_pe {
    input:
      r1_fastq = umiTagger.r1_out_fastq,
      r2_fastq = umiTagger.r3_out_fastq,
      reference = reference
  }

  call sort_and_index {
    input:
      input_bam = hisat2_align_pe.out_bam
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

    File aligned_bam = hisat2_align_pe.out_bam

    File sorted_bam = sort_and_index.output_bam
    File sorted_bam_index = sort_and_index.output_bam_index

    File deduplicated_bam = removeDuplicates.deduplicated_bam
    Array[File] dedup_stats = removeDuplicates.dedup_stats
  }

}


