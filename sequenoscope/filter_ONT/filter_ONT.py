#!/usr/bin/env python
import os
import sys
import time
import argparse as ap
from sequenoscope.utils.__init__ import format_time
from sequenoscope.version import __version__
from sequenoscope.utils.parser import GeneralSeqParser
from sequenoscope.utils.sequence_class import Sequence
from sequenoscope.filter_ONT.seq_summary_processing import SeqSummaryProcesser
from sequenoscope.filter_ONT.seqtk import SeqtkRunner
from sequenoscope.filter_ONT.barcode_statistics import BarcodeStatistics
import warnings
warnings.simplefilter('always', UserWarning)

def parse_args():
    parser = ap.ArgumentParser(prog="sequenoscope",
                               usage="sequenoscope filter_ONT --input_fastq <file.fq> --input_summary <seq_summary.txt> -o <out.fastq> [options]\nFor help use: sequenoscope filter_ONT -h or sequenoscope filter_ONT --help", 
                                description="%(prog)s version {}: a flexible tool for processing multiplatform sequencing data: analyze, subset/filter, compare and visualize.".format(__version__), 
                                formatter_class= ap.RawTextHelpFormatter)

    parser._optionals.title = "Arguments"

    parser.add_argument("--input_fastq", metavar="", nargs="+", help="Path to adaptive sequencing fastq files to process. Not required when using --summarize.")
    parser.add_argument("--input_summary", metavar="", required=True, help="[REQUIRED] Path to ONT sequencing summary file.")
    parser.add_argument("-o", "--output", metavar="", required=True, help="[REQUIRED] Output directory designation")
    parser.add_argument("-op", "--output_prefix", metavar="", default= "sample", help="Output file prefix designation. default is [sample]")
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
    parser.add_argument('--summarize', required=False, help='Generate barcode statistics. Must specify an input summary and output directory', action='store_true')
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
    summarize = args.summarize
    force = args.force

    
    if summarize:
        expected_output_file = os.path.join(out_directory, "{}_barcode_statistics.csv".format(out_prefix))
        if os.path.exists(expected_output_file) and not force:
            print(f"Warning: File {expected_output_file} already exists. Use '--force' option to overwrite.")
            sys.exit()

        if not args.input_summary or not args.output:
            print("Error: When using --summarize, you must specify both --input_summary and --output.")
            sys.exit()
        if args.input_fastq:
            print("-"*40)
            print("WARNING: When using --summarize, the --input_fastq is ignored.")
            print("-"*40)

        print("-"*40)
        print("sequenoscope filter_ONT version {}: Generating barcode statistics...".format(__version__))
        print("-"*40)

        barcode_stats = BarcodeStatistics(input_summary, out_directory, out_prefix)
        barcode_stats.generate_statistics()

        print("Barcode statistics saved to:", barcode_stats.result_files["output_csv_file"])
        print("-"*40)

        sys.exit()
    else:
        if not args.input_fastq:
            print("Error: You must specify --input_fastq for the original workflow.")
            sys.exit()

    params_list = [
        ("Mode", "Filter_ONT"),
        ("Input(s)", ', '.join(input_fastq)),
        ("Outputs folder", out_directory),
        ("Sequenicing summary", input_summary)
    ]

    print("-" * 40)
    print("Input Parameters Summary:")
    print("-" * 40)

    for param_name, param_value in params_list:
        print(f"{param_name}: {param_value}")
       
    print("-"*40)
    print("sequenoscope filter_ONT version {}: Filtering reads...".format(__version__))
    print("-"*40)

    start_time = time.time()

    ## intializing directory for files

    if not os.path.isdir(out_directory):
        os.mkdir(out_directory, 0o755)
    elif not force:
        print("Error directory {} already exists, if you want to overwrite existing results then specify --force".format(out_directory))
        sys.exit()

    ## parsing seq summary file

    print("-"*40)
    print("Extracting reads...")
    print("-"*40)

# Define the required columns based on command-line arguments
    required_columns = [
        "read_id",
        "channel" if min_ch or max_ch else None,
        "start_time" if min_start or max_start else None,
        "duration" if min_dur or max_dur else None,
        "sequence_length_template" if min_len or max_len else None,
        "mean_qscore_template" if min_q or max_q else None,
        "end_reason" if as_class else None,
    ]
    required_columns = [col for col in required_columns if col]  # Remove None values

    # Initialize the parser and parse the summary
    seq_summary_parsed = GeneralSeqParser(input_summary, "seq_summary", required_columns)
    # seq_summary_parsed.required_columns = required_columns
    # seq_summary_parsed.parse_seq_summary()

    ## producing read list

    seq_summary_process = SeqSummaryProcesser(seq_summary_parsed, out_directory, "{}_read_id_list".format(out_prefix), classification= as_class, 
                                            min_ch=min_ch, max_ch=max_ch, min_dur=min_dur, max_dur=max_dur, min_start_time=min_start,
                                            max_start_time=max_start, min_q=min_q, max_q=max_q, min_len=min_len, max_len=max_len)

    
    seq_summary_process.generate_read_ids()

    ## producing fastq via seqtk

    print("-"*40)
    print("Subsetting fastq file...")
    print("-"*40)

    sequencing_sample = Sequence("ONT", input_fastq)
    seqtk_subset = SeqtkRunner(sequencing_sample, seq_summary_process.result_files["filtered_read_id_list"], out_directory, "{}_filtered_fastq".format(out_prefix))
    seqtk_subset.subset_fastq()

    end_time = time.time()
    total_runtime_minutes = end_time - start_time

    print("-"*40)
    print("All Done!")
    print(f"total runtime: {format_time(total_runtime_minutes)}")
    print("-"*40)