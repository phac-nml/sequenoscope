#!/usr/bin/env python
import argparse as ap
import os
import sys
from Sequenoscope.version import __version__
from Sequenoscope.utils.parser import GeneralSeqParser 
from Sequenoscope.utils.sequence_class import Sequence
from Sequenoscope.filter_ONT.seq_summary_processing import SeqSummaryProcesser
from Sequenoscope.filter_ONT.seqtk import SeqtkRunner

def parse_args():
    parser = ap.ArgumentParser(prog="sequenoscope",
                               usage="sequenoscope filter_ONT --input_fastq <file.fq> --input_summary <seq_summary.txt> -o <out.fastq> [options]\nFor help use: sequenoscope filter_ONT -h or sequenoscope filter_ONT --help", 
                                description="%(prog)s version {}: a tool for analyzing and processing sequencing data.".format(__version__), 
                                formatter_class= ap.RawTextHelpFormatter)

    parser._optionals.title = "Arguments"

    parser.add_argument("--input_fastq", metavar="", required=True, nargs="+", help="[REQUIRED] Path to adaptive sequencing fastq files to process.")
    parser.add_argument("--input_summary", metavar="", required=True, help="[REQUIRED] Path to ONT sequencing summary file.")
    parser.add_argument("-o", "--output", metavar="", required=True, help="[REQUIRED] Output directory designation")
    parser.add_argument("-o_pre", "--output_prefix", metavar="", default= "sample", help="Output file prefix designation. default is [sample]")
    parser.add_argument("-cls", "--classification", default= "all", metavar="", type= str, choices=['all', 'unblocked', 'stop_receiving', 'no_decision'], help="a designation of the adaptive-sampling sequencing decision classification ['unblocked', 'stop_receiving', or 'no_decision']")
    parser.add_argument("-min_ch", "--minimum_channel", default= 1, metavar="", type=int, help="a designation of the minimum channel/pore number for filtering reads")
    parser.add_argument("-max_ch", "--maximum_channel", default= 512, metavar="", type=int, help="a designation of the maximum channel/pore number for filtering reads")
    parser.add_argument("-min_dur", "--minimum_duration", default= 0, metavar="", type=float, help="a designation of the minimum duration of the sequencing run in SECONDS for filtering reads")
    parser.add_argument("-max_dur", "--maximum_duration", default= 100, metavar="", type=float, help="a designation of the maximum duration of the sequencing run in SECONDS for filtering reads")
    parser.add_argument("-min_start", "--minimum_start_time", default= 0, metavar="", type=float, help="a designation of the minimum start time of the sequencing run in SECONDS for filtering reads")
    parser.add_argument("-max_start", "--maximum_start_time", default= 259200,metavar="", type=float, help="a designation of the maximum start time of the sequencing run in SECONDS for filtering reads")
    parser.add_argument("-min_q", "--minimum_q_score", metavar="", default= 0, type=int, help="a designation of the minimum q score for filtering reads")
    parser.add_argument("-max_q", "--maximum_q_score", metavar="", default= 100, type=int, help="a designation of the maximum q score for filtering reads")
    parser.add_argument("-min_len", "--minimum_length", metavar="", default= 0, type=int, help="a designation of the minimum read length for filtering reads")
    parser.add_argument("-max_len", "--maximum_length", metavar="", default= 50000,type=int, help="a designation of the maximum read length for filtering reads")
    parser.add_argument('--force', required=False, help='Force overwite of existing results directory', action='store_true')
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + __version__)
    return parser.parse_args()

def run():
    args = parse_args()
    input_fastq = args.input_fastq
    input_summary = args.input_summary
    out_directory = args.output
    out_prefix = args.output_prefix
    as_class = args.classification
    min_ch = args.minimum_channel
    max_ch = args.maximum_channel
    min_dur = args.minimum_duration
    max_dur = args.maximum_duration
    min_start = args.minimum_start_time
    max_start = args.maximum_start_time
    min_q = args.minimum_q_score
    max_q = args.maximum_q_score
    min_len = args.minimum_length
    max_len = args.maximum_length
    force = args.force

    print("-"*40)
    print("Sequenoscope filter_ONT version {}: filtering ONT reads based on paramters".format(__version__))
    print("-"*40)

    ## intializing directory for files

    if not os.path.isdir(out_directory):
        os.mkdir(out_directory, 0o755)
    elif not force:
        print("Error directory {} already exists, if you want to overwrite existing results then specify --force".format(out_directory))
        sys.exit()

    ## parsing seq summary file

    print("-"*40)
    print("Processing seq summary file and extracting reads based on input parameters...")
    print("-"*40)

    seq_summary_parsed = GeneralSeqParser(input_summary, "seq_summary")

    ## producing read list

    seq_summary_process = SeqSummaryProcesser(seq_summary_parsed, out_directory, "{}_read_id_list".format(out_prefix), classification= as_class, 
                                            min_ch=min_ch, max_ch=max_ch, min_dur=min_dur, max_dur=max_dur, min_start_time=min_start,
                                            max_start_time=max_start, min_q=min_q, max_q=max_q, min_len=min_len, max_len=max_len)

    
    seq_summary_process.generate_read_ids()

    ## producing fastq via seqtk

    print("-"*40)
    print("Subsetting fastq file based on extracted reads...")
    print("-"*40)

    sequencing_sample = Sequence("ONT", input_fastq)
    seqtk_subset = SeqtkRunner(sequencing_sample, seq_summary_process.result_files["filtered_read_id_list"], out_directory, "{}_filtered_fastq".format(out_prefix))
    seqtk_subset.subset_fastq()

    print("-"*40)
    print("All Done!")
    print("-"*40)