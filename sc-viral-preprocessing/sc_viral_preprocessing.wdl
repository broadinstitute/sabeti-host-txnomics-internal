version 1.0

task SortAndIndex {
  input {
    File input_bam
    String docker  = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "sorted.bam"
  String output_bai_name = "sorted.bai"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 3 + 10)
  Int preemptible =3 

  command {
    set -e

    java -Xmx6G -jar $PICARD_JAR_PATH SortSam I=~{input_bam} O=~{output_bam_name} SORT_ORDER=coordinate
    java -jar $PICARD_JAR_PATH BuildBamIndex I=~{output_bam_name}
  }

  runtime {
    docker: docker
    memory: "8 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = output_bam_name
    File output_bai = output_bai_name
  }

}

task FilterReads {
  input {
    File input_bam
    File input_bam_index
    String region_keep
    
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "filtered.bam"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e
    
    samtools view -b ~{input_bam} ~{region_keep} > ~{output_bam_name}
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
    File output_bam = output_bam_name
  }

}

task SamToFastq {
  input {
    File input_bam
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_fastq_name = "reverted.fastq"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 4 + 10)
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
    File output_fastq = output_fastq_name
  }
}

task BWAalign {
  input {
    File input_fastq
    File reference_bundle
    String reference_bundle_prefix
    
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-bwa:0.7.17"
  }

  String output_bam_name = "aligned.bam"
  Int cpu = 1
  Int disk = ceil(size(input_fastq, "GiB") * 4 + 10)
  Int preemptible = 3

  command {
    set -e

    tar xvzf ~{reference_bundle}
    
    bwa mem -t ~{cpu} ~{reference_bundle_prefix} ~{input_fastq} | samtools view -b - > ~{output_bam_name}
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

task MergeBamAlignment {
  input {
    File input_aligned_bam
    File input_unaligned_bam
    File reference_fasta
    
    String docker="us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "merged_alignments.bam"
  Int cpu = 1
  Int disk = ceil(size(input_aligned_bam, "GiB") * 4 + size(input_unaligned_bam, "GiB") + 10)
  Int preemptible = 3

  command {
    set -e
    java -jar $PICARD_JAR_PATH CreateSequenceDictionary R=~{reference_fasta}
    java -jar $PICARD_JAR_PATH MergeBamAlignment ALIGNED=~{input_aligned_bam} UNMAPPED=~{input_unaligned_bam} O=~{output_bam_name} REFERENCE_SEQUENCE=~{reference_fasta}
  }

  output {
    File output_bam = output_bam_name
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disk: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
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
    Int disk = ceil((size(bam_input, "GiB") + size(annotations_gtf, "GiB")) * 3) + 20
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

task FilterUnmapped {
  input {
    File bam_input
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"
  }

  String output_bam_name = "output.bam"
  Int cpu = 1
  Int disk = ceil(size(bam_input, "GiB")  * 4 + 10)
  Int preemptible = 3

  command {
    set -e

    samtools view -b -F 4 ~{bam_input} > ~{output_bam_name}
  }

  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = output_bam_name
  }
}

task umiCollapser {
  input {
    File bam_input
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-umi_collapser:0.0.1"
  }

  String output_bam_name = "output.bam"
  Int machine_mem_mb = 8250
  Int cpu = 1
  Int disk = ceil(size(bam_input, "GiB") * 3 + 10)
  Int preemptible = 3

  command {
    set -e
    
    umi_collapser -i ~{bam_input} -o ~{output_bam_name} --cell_barcode_tag XC --molecular_barcode_tag XM --gene_tag gn --calling_method posterior --verbose
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = output_bam_name
  }
}

workflow scViralPreprocess {
  input {
    File input_bam
    String viral_contig
    File reference_bundle
    String reference_bundle_prefix
    File annotations_gtf
    File reference_fasta
  }

  String version = "scViralPreprocess_v0.0.1"

  call SortAndIndex {
    input:
      input_bam = input_bam
  }

  call FilterReads {
    input:
      input_bam = SortAndIndex.output_bam,
      input_bam_index = SortAndIndex.output_bai,
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
      input_fastq = SamToFastq.output_fastq,
      reference_bundle = reference_bundle,
      reference_bundle_prefix = reference_bundle_prefix
  }

  call MergeBamAlignment {
    input:
      input_aligned_bam = BWAalign.output_bam,
      input_unaligned_bam = RevertSam.output_bam,
      reference_fasta = reference_fasta
  }

  call TagReadWithGeneFunction {
    input:
      bam_input = MergeBamAlignment.output_bam,
      annotations_gtf = annotations_gtf
  }

  call FilterUnmapped {
    input:
      bam_input = TagReadWithGeneFunction.bam_output
  } 

  call umiCollapser {
    input:
      bam_input = FilterUnmapped.bam_output
  }

  call SortAndIndex as FinalSortAndIndex {
    input:
      input_bam = umiCollapser.bam_output
  }

  output {
    String pipeline_version = version
    File output_bam = FinalSortAndIndex.output_bam
  }
}