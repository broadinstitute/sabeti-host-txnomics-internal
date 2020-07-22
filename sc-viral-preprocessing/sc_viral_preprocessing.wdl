version 1.0

task FilterReads {
  command {
    samtools view ~{input_bam} ~{region_keep}
  }
}

task RevertSam {
  command {
    java -jar $PICARD_JAR_PATH RevertSam I=${INPUT_BAM} O=reverted.bam
  }
}

task SamToFastq {
  command {
    java -jar $PICARD_JAR_PATH SamToFastq I=reverted.bam FASTQ=reverted.fastq
  }
}

task BWAalign {
  command {
    bwa mem -t 12 bwa_ebov/kitwit9510621_KU182905.1.fasta reverted.fastq | samtools view -b - > realigned.bam
  }

}

task MergeBamAlignment {
  command {
    java -jar $PICARD_JAR_PATH MergeBamAlignment ALIGNED=realigned.bam UNMAPPED=reverted.bam O=merge_alignments.bam REFERENCE_SEQUENCE=bwa_ebov/kitwit9510621_KU182905.1.fasta
  }
}

task TagReadWithGeneFunction {
  input {
    File annotations_gtf
    File bam_input

    String gene_name_tag = "gn"
    String gene_strand_tag = "gs"
    String gene_function_tag = "gf"

    String use_strand_info = "true"

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:2.3.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil((size(bam_input, "Gi") + size(annotations_gtf, "Gi")) * 3) + 20
    Int preemptible = 3
  }

  meta {
    description: "Tags any read in bam_input that overlaps an intron or exon interval with the gene that those interals correspond to."
  }

  parameter_meta {
    annotations_gtf: "GTF annotation file for the species that bam input is derived from. Each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id, no white space at the end of any record and must be in gtf format."
    bam_input: "Aligned bam file."
    gene_name_tag: "the tag used to denote gene name in the bam (default: gn)"
    gene_strand_tag: "the tag used to denote gene strand in the bam (default: gs)"
    gene_function_tag: "the tag used to denote gene function (INTRONIC, EXONIC, ...) in the output bam (default: gf)"

    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

 command {
    set -e

    TagReadWithGeneFunction \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      GENE_NAME_TAG=${gene_name_tag} \
      GENE_STRAND_TAG=${gene_strand_tag} \
      GENE_FUNCTION_TAG=${gene_function_tag} \
      SUMMARY=gene_exon_tag_summary.log \
      ANNOTATIONS_FILE=${annotations_gtf} \
      USE_STRAND_INFO=${use_strand_info}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 


orkflow scViralPreprocess {
  input {
    File input_bam
    String viral_contig
    File reference
  }

  String version = "scViralPreprocess_v0.0.1"

  call FilterReads {
    input:
      input_bam = input_bam
      viral_contig = viral_contig
  }

  call RevertSam {
    input:
      input_bam = FilterReads.output_bam
  }

  call SamToFastq {
    input:
      input_bam = RevertSam.output_bam
  }

  call BWAalign {
    input:
      input_bam = SamToFastq.output_bam
      reference = reference
  }

  call MergeBamAlignment {
    input:
      input_aligned = BWAalign.output_bam
      input_unmapped = FilterReads.output_bam
  }

  call TagReadWithGeneFunction {
    input:
      bam_input = MergeBamAlignment.output_bam
  }
}