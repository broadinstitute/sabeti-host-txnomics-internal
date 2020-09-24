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

task StringTie2 {
  input {
    File input_bam
    File annotation_gtf
  }

  String output_annotation_filename = "output_annotation.gtf"
  String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-stringtie2:0.0.1"
  Int cpu = 1
  Int disk = ceil(size(input_bam, "GiB") * 2 + 10)
  Int preemptible = 3

  command {
    set -e
    stringtie -G ~{annotation_gtf} -p ~{cpu} -o ~{output_annotation_filename} ~{input_bam}
  }

  runtime {
    docker: docker
    memory: "8 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_annotation = output_annotation_filename
  }
}

task StringTie2_merge {
  input {
    Array[File] input_gtf
    File guide_gtf
  }

  String output_annotation_filename = "merged_annot.gtf"
  String docker = "us.gcr.io/broad-dsde-methods/sabeti-bulk-plp-stringtie2:0.0.1"
  Int cpu = 1
  Int disk = 500
  Int preemptible = 3

  command {
    set -e
    stringtie --merge -G ~{guide_gtf} -o ~{output_annotation_filename} -m 50 -c 0 -F 0 -T 0 -f 0.01 -l MSTRG ~{sep=' ' input_gtf}
  }

  runtime {
    docker: docker
    memory: "8 GiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File merged_annotation = output_annotation_filename
  }
} 

workflow buildAnnotationNoMerge {
  input {
    File input_bam
    File annotations_gtf
  }

  String version = "build-annotation_v0.0.1"

  call SortAndIndex as ScatterSort {
    input:
      input_bam = input_bam
  }

  call StringTie2 {
    input:
      input_bam = ScatterSort.output_bam,
      annotation_gtf = annotations_gtf
  }

  output {
     String pipeline_version = version
     File stringtie_output_gtf = StringTie2.output_annotation
  }
} 