#!/usr/bin/env python
import os
import sys
import time
import logging
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
    parser = ap.ArgumentParser(
        prog="sequenoscope",
        usage="sequenoscope filter_ONT --input_fastq <file.fq> --input_summary <seq_summary.txt> -o <out.fastq> [options]\nFor help use: sequenoscope filter_ONT -h or sequenoscope filter_ONT --help", 
        description="%(prog)s version {}: a flexible tool for processing multiplatform sequencing data: analyze, subset/filter, compare and visualize.".format(__version__), 
        formatter_class= ap.RawTextHelpFormatter
    )

    parser._optionals.title = "Arguments"

    parser.add_argument("--input_fastq", metavar="", nargs="+", help="Path to adaptive sequencing fastq files to process. Not required when using --summarize.")
    parser.add_argument("--input_summary", metavar="", required=True, help="[REQUIRED] Path to ONT sequencing summary file.")
    parser.add_argument("-o", "--output", metavar="", required=True, help="[REQUIRED] Output directory designation")
    parser.add_argument("-op", "--output_prefix", metavar="", default="sample", help="Output file prefix designation. default is [sample]")
    parser.add_argument("-cls", "--classification", default="all", metavar="", type=str, 
                        choices=['all', 'unblocked', 'stop_receiving', 'no_decision'],
                        help=       "Adaptive-sampling sequencing decision classification based on the 'end_reason' column "
                                    "in the sequencing summary file. Options:\n"
                                    "  - 'unblocked': 'data_service_unblock_mux_change'\n"
                                    "  - 'stop_receiving': 'signal_positive'\n"
                                    "  - 'no_decision': ('signal_negative', 'unblock_mux_change')\n"
                                    "  - 'all': Includes all classifications.\n"
                                    "Default: 'all'.")
    parser.add_argument("-min_ch", "--minimum_channel", default=1, metavar="", type=int, help="Minimum channel number for filtering reads. Default=1.")
    parser.add_argument("-max_ch", "--maximum_channel", default=512, metavar="", type=int, help="Maximum channel number for filtering reads. Default=512.")
    parser.add_argument("-min_dur", "--minimum_duration", default=0, metavar="", type=float, help="Minimum duration (s) for filtering reads. Default=0.")
    parser.add_argument("-max_dur", "--maximum_duration", default=100, metavar="", type=float, help="Maximum duration (s) for filtering reads. Default=100.")
    parser.add_argument("-min_start", "--minimum_start_time", default=0, metavar="", type=float, help="Minimum start time (s) for filtering reads. Default=0.")
    parser.add_argument("-max_start", "--maximum_start_time", default=259200, metavar="", type=float, help="Maximum start time (s) for filtering reads. Default=259200.")
    parser.add_argument("-min_q", "--minimum_q_score", metavar="", default=0, type=int, help="Minimum Q score for filtering reads. Default=0.")
    parser.add_argument("-max_q", "--maximum_q_score", metavar="", default=100, type=int, help="Maximum Q score for filtering reads. Default=100.")
    parser.add_argument("-min_len", "--minimum_length", metavar="", default=0, type=int, help="Minimum read length for filtering. Default=0.")
    parser.add_argument("-max_len", "--maximum_length", metavar="", default=50000, type=int, help="Maximum read length for filtering. Default=50000.")
    parser.add_argument('--force', required=False, help='Force overwrite of existing results directory', action='store_true')
    parser.add_argument('--summarize', required=False, action='store_true', help=   "Generate barcode statistics only. This mode works exclusively with the\n"
                                                                                    "'--input_summary' argument. You must specify both '--input_summary' and\n"
                                                                                    "'--output' when using this option. The output is a CSV file containing "
                                                                                    "barcode-level statistics, including:\n"
                                                                                    "  - Barcode ID\n"
                                                                                    "  - Total number of reads per barcode\n"
                                                                                    "  - Minimum, maximum, and mean channel numbers\n"
                                                                                    "  - Minimum, maximum, and mean start times (s)\n"
                                                                                    "  - Minimum, maximum, and mean read durations (s)\n"
                                                                                    "  - Minimum, maximum, and mean sequence lengths (bp)\n"
                                                                                    "  - Minimum, maximum, and mean quality scores (Q).\n"
                                                                                    "This file provides a detailed summary of sequencing metrics per barcode.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + __version__)
    return parser.parse_args()

def count_fastq_reads(fastq_files):
    """Count total number of reads in one or more FASTQ files."""
    total_reads = 0
    for fq in fastq_files:
        with open(fq, 'r') as f:
            # Each read occupies 4 lines in a standard FASTQ file
            lines = sum(1 for _ in f)
            total_reads += lines // 4
    return total_reads

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

    # Set up output directory and logging
    if not os.path.isdir(out_directory):
        os.mkdir(out_directory, 0o755)
    elif not force:
        print(f"Error directory {out_directory} already exists, if you want to overwrite existing results then specify --force", file=sys.stderr)
        sys.exit()

    log_filename = "filter.log"
    log_filepath = os.path.join(out_directory, log_filename)

    logger = logging.getLogger("sequenoscope_filter")
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_filepath, mode='w')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("Starting 'sequenoscope filter_ONT' module.")
    logger.info(f"Version: {__version__}")
    logger.info("Validating input parameters...")

    if summarize:
        logger.info("Summarize mode selected. Generating barcode statistics only.")
        expected_output_file = os.path.join(out_directory, f"{out_prefix}_barcode_statistics.csv")

        if os.path.exists(expected_output_file) and not force:
            logger.error(f"File {expected_output_file} already exists and overwrite not forced.")
            print(f"Warning: File {expected_output_file} already exists. Use '--force' option to overwrite.")
            sys.exit()

        if not args.input_summary or not args.output:
            logger.error("Missing required input_summary or output for summarize mode.")
            print("Error: When using --summarize, you must specify both --input_summary and --output.")
            sys.exit()

        if args.input_fastq:
            logger.info("Input_fastq provided but will be ignored since summarize mode is active.")

        print("-"*40)
        print(f"sequenoscope filter_ONT version {__version__}: Generating barcode statistics...")
        print("-"*40)

        logger.info("Generating barcode statistics...")
        barcode_stats = BarcodeStatistics(input_summary, out_directory, out_prefix)
        barcode_stats.generate_statistics()
        logger.info(f"Barcode statistics saved to {barcode_stats.result_files['output_csv_file']}")

        print("Barcode statistics saved to:", barcode_stats.result_files["output_csv_file"])
        print("-"*40)

        logger.info("Summarize mode completed successfully. Exiting.")
        sys.exit()
    else:
        if not input_fastq:
            logger.error("input_fastq is required for non-summarize mode.")
            print("Error: You must specify --input_fastq for the original workflow.")
            sys.exit()

    # Count input reads
    input_read_count = count_fastq_reads(input_fastq)

    # Log input parameters
    logger.info("Input Parameters:")
    logger.info("-" * 40)
    logger.info("Mode: Filter_ONT")
    logger.info(f"Input FASTQ(s): {', '.join(input_fastq)}")
    logger.info(f"Output directory: {out_directory}")
    logger.info(f"Sequencing summary: {input_summary}")
    logger.info(f"Classification filter: {as_class}")
    logger.info(f"Min channel: {min_ch}, Max channel: {max_ch}")
    logger.info(f"Min duration: {min_dur}, Max duration: {max_dur}")
    logger.info(f"Min start time: {min_start}, Max start time: {max_start}")
    logger.info(f"Min Q score: {min_q}, Max Q score: {max_q}")
    logger.info(f"Min length: {min_len}, Max length: {max_len}")
    logger.info(f"Total input reads (approx.): {input_read_count}")
    logger.info("-" * 40)
    logger.info("All input parameters validated successfully.")

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
    print(f"sequenoscope filter_ONT version {__version__}: Filtering reads...")
    print("-"*40)

    start_time = time.time()

    logger.info("Filtering reads based on specified parameters.")

    print("-"*40)
    print("Extracting reads...")
    print("-"*40)
    logger.info("Extracting required columns from sequencing summary and filtering read IDs.")

    required_columns = [
        "read_id",
        "channel" if (min_ch or max_ch) else None,
        "start_time" if (min_start or max_start) else None,
        "duration" if (min_dur or max_dur) else None,
        "sequence_length_template" if (min_len or max_len) else None,
        "mean_qscore_template" if (min_q or max_q) else None,
        "end_reason" if as_class else None,
    ]
    required_columns = [col for col in required_columns if col]

    seq_summary_parsed = GeneralSeqParser(input_summary, "seq_summary", required_columns)

    seq_summary_process = SeqSummaryProcesser(
        seq_summary_parsed, out_directory, f"{out_prefix}_read_id_list", 
        classification=as_class, 
        min_ch=min_ch, max_ch=max_ch, 
        min_dur=min_dur, max_dur=max_dur, 
        min_start_time=min_start, max_start_time=max_start, 
        min_q=min_q, max_q=max_q, 
        min_len=min_len, max_len=max_len
    )

    seq_summary_process.generate_read_ids()
    logger.info("Filtered read ID list generated successfully.")

    # Count how many reads passed filters
    filtered_read_id_list_file = seq_summary_process.result_files["filtered_read_id_list"]
    with open(filtered_read_id_list_file, 'r') as f:
        filtered_read_count = sum(1 for _ in f)

    print("-"*40)
    print("Subsetting fastq file...")
    print("-"*40)
    logger.info("Subsetting FASTQ file using seqtk.")

    sequencing_sample = Sequence("ONT", input_fastq)
    seqtk_subset = SeqtkRunner(sequencing_sample, filtered_read_id_list_file, out_directory, f"{out_prefix}_filtered_fastq")
    seqtk_subset.subset_fastq()
    logger.info("FASTQ subsetting completed successfully.")

    end_time = time.time()
    total_runtime_minutes = end_time - start_time

    # Print and log summary of filtering
    reads_removed = input_read_count - filtered_read_count
    logger.info(f"Input reads: {input_read_count}")
    logger.info(f"Filtered reads: {filtered_read_count}")
    logger.info(f"Reads removed: {reads_removed}")

    print("-"*40)
    print("All Done!")
    print(f"Total input reads: {input_read_count}")
    print(f"Filtered reads (passing criteria): {filtered_read_count}")
    print(f"Reads removed: {reads_removed}")
    print(f"total runtime: {format_time(total_runtime_minutes)}")
    print("-"*40)

    logger.info("Filtering pipeline completed successfully.")
    logger.info(f"Total runtime: {format_time(total_runtime_minutes)}")
    logger.info("All operations are complete.")
