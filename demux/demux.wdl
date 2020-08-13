version 1.0

task demux {
    input {
        File flowcell
        String flowcell_basecall_path
        File sample_sheet
        String flowcell_barcode
        String machine_name
        String run_barcode
        String lane
        String read_structure
    }

    Int cpu = 1
    Int disk = ceil(size(flowcell,"GiB") * 10 + 10)
    Int preemptible = 3
    String docker = "us.gcr.io/broad-dsde-methods/sabeti-picard:2.23.2"

    command {
        set -e

        mkdir flowcell_dir

        tar -C flowcell_dir -xzf ~{flowcell}

        printf "library_name\tbarcode_sequence_1\tbarcode_sequence_2\n" > barcodes.txt
        cat ~{sample_sheet} >> barcodes.txt

        java -jar $PICARD_JAR_PATH ExtractIlluminaBarcodes \
            BASECALLS_DIR=flowcell_dir/~{flowcell_basecall_path} \
            LANE=~{lane} \
            READ_STRUCTURE=~{read_structure} \
            BARCODE_FILE=barcodes.txt \
            METRICS_FILE=metrics_output.txt

        printf "OUTPUT_PREFIX\tBARCODE_1\tBARCODE_2\n" > multiplex_params.txt
        cat ~{sample_sheet} >> multiplex_params.txt

        mkdir output
        cd output || exit

        java -jar $PICARD_JAR_PATH IlluminaBasecallsToFastq \
            READ_STRUCTURE=~{read_structure} \
            BASECALLS_DIR=../flowcell_dir/~{flowcell_basecall_path} \
            LANE=~{lane} \
            FLOWCELL_BARCODE=~{flowcell_barcode} \
            MACHINE_NAME=~{machine_name} \
            RUN_BARCODE=~{run_barcode} \
            MULTIPLEX_PARAMS=../multiplex_params.txt \
            INGORE_UNEXPECTED=true \
            COMPRESS_OUTPUTS=true

    }

    runtime {
        docker: docker
        disks: "local-disk ~{disk} HDD"
        memory: "8 GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        Array[File] fastq_out = glob("output/*.fastq.gz")
        File barcode_metrics ="metrics_output.txt"
    }
}

workflow demultiplex {
    input {
        File flowcell
        File sample_sheet
        String flowcell_basecall_path = ""
        String flowcell_barcode
        String machine_name
        String run_barcode
        String lane
        String read_structure
    }

    String version = "demux_v0.0.1"

    call demux {
        input:
            flowcell = flowcell,
            sample_sheet = sample_sheet,
            flowcell_basecall_path = flowcell_basecall_path,
            flowcell_barcode = flowcell_barcode,
            machine_name = machine_name,
            run_barcode = run_barcode,
            lane = lane,
            read_structure = read_structure
    }

    output {
        Array[File] fastq_out = demux.fastq_out
        String pipeline_version = version
    }

}
