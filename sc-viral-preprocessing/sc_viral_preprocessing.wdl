version 1.0

task FilterReads {
  input {
    File input_bam
    String region_keep
    
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "filtered.bam"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e
    
    samtools view ~{input_bam} ~{region_keep} > ~{output_bam_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = output_bam_name
  }
}

task RevertSam {
  input {
    File input_bam
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "reverted.bam"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e
    
    java -jar $PICARD_JAR_PATH RevertSam I=~{input_bam} O=~{output_bam_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = ~{output_bam_name}
  }

}

task SamToFastq {
  input {
    File input_bam
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_fastq_name = "reverted.fastq"
  Int cpu = 1
  Int disk = ceil(size(input_bam), "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e 

    java -jar $PICARD_JAR_PATH SamToFastq I=~{input_bam} FASTQ=~{output_fastq_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    output_fastq = ~{output_fastq_name}
  }
}

task BWAalign {
  input {
    File input_fastq
    File reference_bundle
    
    String docker = ""
  }

  String output_bam_name = "aligned.bam"
  Int cpu = 1
  Int disk = ceil(size(input_fastq, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e
    
    bwa mem -t ~{cpu} bwa_ebov/kitwit9510621_KU182905.1.fasta ~{input_fastq} | samtools view -b - > ~{output_bam_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = ~{output_bam_name}
  }
}

task MergeBamAlignment {
  input {
    File input_aligned_bam
    File input_unaligned_bam
    File reference_fasta
    String docker="us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "merged_alignments.bam"

  command {
    java -jar ~{PICARD_JAR_PATH} MergeBamAlignment ALIGNED=~{input_aligned_bam} UNMAPPED=~{input_unaligned_bam} O=~{output_bam_name} REFERENCE_SEQUENCE=~{reference_fasta}
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


workflow scViralPreprocess {
  input {
    File input_bam
    String viral_contig
    File reference_bundle
    File annotations_gtf
  }

  String version = "scViralPreprocess_v0.0.1"

  call FilterReads {
    input:
      input_bam = input_bam
      region_keep = viral_contig
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
      reference_bundle = reference_bundle
  }

  call MergeBamAlignment {
    input:
      input_aligned = BWAalign.output_bam
      input_unmapped = RevertSam.output_bam
  }

  call TagReadWithGeneFunction {
    input:
      bam_input = MergeBamAlignment.output_bam
      annotations_gtf = annotations_gtf
  }

  output {
    String pipeline_version = version
    File output_bam = TagReadWithGeneFunction
  }
}