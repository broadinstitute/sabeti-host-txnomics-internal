#!/usr/bin/env python3

import gzip
from Bio import SeqIO
import argparse



def proc_file(r1_filename, r2_filename, r3_filename, r1_out_filename, r3_out_filename):
    i = 0

    with gzip.open(r1_filename,'rt') as r1_in:
        with gzip.open(r2_filename, 'rt') as r2_in:
            with gzip.open(r3_filename, 'rt') as r3_in:
                with gzip.open(r1_out_filename, 'wt') as r1_out:
                    with gzip.open(r3_out_filename, 'wt') as r3_out:
                        r1_parser = SeqIO.parse(r1_in, 'fastq')
                        r2_parser = SeqIO.parse(r2_in, 'fastq')
                        r3_parser = SeqIO.parse(r3_in, 'fastq')

                        for rec1 in r1_parser:
                            rec2 = r2_parser.__next__()
                            rec3 = r3_parser.__next__()

                            umi = rec2.seq

                            rec1.id = rec1.id + '_' + umi.__str__()
                            rec3.id = rec3.id + '_' + umi.__str__()

                            SeqIO.write(rec1, r1_out, 'fastq')
                            SeqIO.write(rec3, r3_out, 'fastq')

                            i += 1
                            if (i % 10000 == 0):
                                print(i)

def main():
    description = "Append R2 sequence as UMI tag to read ids"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--r1',
                        required=True,
                        help="read 1 input FASTQ")
    parser.add_argument('--r2',
                        required=True,
                        help="read 2 input FASTQ")
    parser.add_argument('--r3',
                        required=True,
                        help="read 3 input FASTQ")
    parser.add_argument('--r1-output',
                        required=True,
                        help="read 1 output FASTQ")
    parser.add_argument('--r3-output',
                        required=True,
                        help="read 3 output FASTQ")

    args = parser.parse_args()


    proc_file(args.r1, args.r2, args.r3, args.r1_output, args.r3_output)

    return 0


if __name__ == "__main__":
    main()
