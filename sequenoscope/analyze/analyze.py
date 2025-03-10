#!/usr/bin/env python
import os
import sys
import time
import logging
import argparse as ap
import gzip
import shutil

from sequenoscope.utils.__init__ import format_time
from sequenoscope.constant import SequenceTypes
from sequenoscope.version import __version__
from sequenoscope.utils.parser import GeneralSeqParser 
from sequenoscope.utils.sequence_class import Sequence
from sequenoscope.analyze.minimap2 import Minimap2Runner
from sequenoscope.analyze.fastP import FastPRunner
from sequenoscope.analyze.processing import SamBamProcessor
from sequenoscope.analyze.fastq_extractor import FastqExtractor
from sequenoscope.analyze.seq_manifest import SeqManifest
from sequenoscope.utils.parser import FastqPairedEndRenamer
from sequenoscope.analyze.seq_manifest import SeqManifestSummary
from sequenoscope.analyze.mash import MashSketcher
import warnings
warnings.simplefilter('always', UserWarning)


def decompress_gz_file(gz_path, dest_dir):
    """
    Decompress a .gz file into dest_dir.
    
    Returns the path to the decompressed file.
    """
    base_name = os.path.basename(gz_path)
    if base_name.endswith('.gz'):
        out_file = os.path.join(dest_dir, base_name[:-3])
        with gzip.open(gz_path, 'rb') as f_in, open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return out_file
    return gz_path


def parse_args():
    parser = ap.ArgumentParser(
        prog="sequenoscope",
        usage="sequenoscope analyze --input_fastq <file.fq> --input_reference <ref.fasta> -o <out> -seq_type <sr> [options]\nFor help use: sequenoscope analyze -h or --help",
        description="%(prog)s version {}: a flexible tool for processing multiplatform sequencing data: analyze, subset/filter, compare and visualize.".format(__version__),
        formatter_class=ap.RawTextHelpFormatter
    )

    # USER OPTIONS: Essential inputs and outputs
    user_group = parser.add_argument_group("USER OPTIONS", "Direct input files and basic parameters.")
    user_group.add_argument("--input_fastq", metavar="", required=True, nargs="+",
                        help="[REQUIRED] Path to 1 (SE) or 2 (PE) FASTQ files to process.")
    user_group.add_argument("--input_reference", metavar="", required=True,
                        help="[REQUIRED] Path to a single reference FASTA file.")
    user_group.add_argument("-seq_sum", "--sequencing_summary", metavar="",
                        help="(Optional) Path to sequencing summary for manifest creation.")
    user_group.add_argument("-o", "--output", metavar="", required=True,
                        help="[REQUIRED] Output directory designation.")
    user_group.add_argument("-op", "--output_prefix", metavar="", default="sample",
                        help="Output file prefix designation. Default is 'sample'.")
    user_group.add_argument("-seq_type", "--sequencing_type", required=True, metavar="", type=str, choices=['SE', 'PE'],
                        help="[REQUIRED] Sequencing type: 'SE' for single-end or 'PE' for paired-end.")

    # FILTER OPTIONS: Parameters to filter/trim FASTQ reads.
    filter_group = parser.add_argument_group("FILTER OPTIONS", "Parameters to filter/trim FASTQ reads.")
    filter_group.add_argument("-min_cov", "--minimum_coverage", default=1, metavar="", type=int,
                        help="Minimum coverage threshold; default is 1.")
    filter_group.add_argument("-t", "--threads", default=1, metavar="", type=int,
                        help="Number of threads to use.")
    filter_group.add_argument("-min_len", "--minimum_read_length", default=15, metavar="", type=int,
                        help="Minimum read length; default is 15.")
    filter_group.add_argument("-max_len", "--maximum_read_length", default=0, metavar="", type=int,
                        help="Maximum read length; default is 0 (no limit).")
    filter_group.add_argument("-trm_fr", "--trim_front_bp", default=0, metavar="", type=int,
                        help="Bases to trim from the front; default is 0.")
    filter_group.add_argument("-trm_tail", "--trim_tail_bp", default=0, metavar="", type=int,
                        help="Bases to trim from the tail; default is 0.")
    filter_group.add_argument("-q", "--quality_threshold", default=15, metavar="", type=int,
                        help="Quality score threshold; default is 15.")

    # Note: Start and end times are fixed internally to 0 and 100 when no sequencing summary is provided.
    # Note: The minimap2 kmer option has been removed; kmer size defaults to 15.
    parser.add_argument('--force', action='store_true', help="Force overwrite of existing results directory.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + __version__)
    return parser.parse_args()


def run():
    args = parse_args()
    input_fastq = args.input_fastq
    input_reference = args.input_reference
    seq_summary = args.sequencing_summary
    out_directory = args.output
    out_prefix = args.output_prefix
    seq_class = args.sequencing_type
    threads = args.threads
    minimap_kmer_size = 15  # default value; option removed
    min_len = args.minimum_read_length
    max_len = args.maximum_read_length
    trim_front = args.trim_front_bp
    trim_tail = args.trim_tail_bp
    quality_threshold = args.quality_threshold
    min_cov = args.minimum_coverage
    force = args.force

    # Fixed default times when no sequencing summary is provided.
    start_time_default = 0
    end_time_default = 100

    # Setup output directory (final outputs: manifests and log remain in out_directory)
    if not os.path.isdir(out_directory):
        os.mkdir(out_directory, 0o755)
    elif not force:
        print(f"Error: Directory {out_directory} already exists. Use --force to overwrite.", file=sys.stderr)
        sys.exit()

    # Create an intermediates subdirectory for all intermediary files.
    intermediate_dir = os.path.join(out_directory, "intermediates")
    if not os.path.isdir(intermediate_dir):
        os.mkdir(intermediate_dir, 0o755)

    # Check and decompress gz FASTQ files if necessary.
    decompressed_fastq = []
    for fq in input_fastq:
        if fq.endswith(".gz"):
            # Decompress into intermediate folder.
            decompressed = decompress_gz_file(fq, intermediate_dir)
            decompressed_fastq.append(decompressed)
        else:
            decompressed_fastq.append(fq)
    input_fastq = decompressed_fastq

    # Setup logging (log file remains in out_directory)
    log_filepath = os.path.join(out_directory, "analyze.log")
    logger = logging.getLogger("sequenoscope_analyze")
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_filepath, mode='w')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("Starting 'sequenoscope analyze' module.")
    logger.info(f"Version: {__version__}")
    logger.info("Validating input parameters...")

    # Validate FASTQ count based on sequencing type.
    if seq_class.upper() == SequenceTypes.paired_end and len(input_fastq) != 2:
        logger.error("Paired-end specified but did not receive exactly 2 FASTQ files.")
        print("Error: Paired-end sequencing requires exactly 2 FASTQ files.", file=sys.stderr)
        sys.exit()
    if seq_class.upper() == SequenceTypes.single_end and len(input_fastq) != 1:
        logger.error("Single-end specified but received multiple FASTQ files.")
        print("Error: Single-end sequencing requires exactly 1 FASTQ file.", file=sys.stderr)
        sys.exit()

    logger.info("Input Parameters:")
    logger.info("-" * 40)
    logger.info(f"Mode: analyze")
    logger.info(f"Input FASTQ: {', '.join(input_fastq)}")
    logger.info(f"Reference FASTA: {input_reference}")
    logger.info(f"Output directory: {out_directory}")
    logger.info(f"Output prefix: {out_prefix}")
    logger.info(f"Sequencing Type: {seq_class}")
    if seq_summary:
        logger.info(f"Sequencing Summary: {seq_summary}")
    logger.info(f"Threads: {threads}")
    logger.info(f"Minimum read length: {min_len}")
    logger.info(f"Maximum read length: {max_len}")
    logger.info(f"Trim front bases: {trim_front}")
    logger.info(f"Trim tail bases: {trim_tail}")
    logger.info(f"Quality threshold: {quality_threshold}")
    logger.info(f"Minimum coverage: {min_cov}")
    logger.info(f"Minimap2 kmer size (default): {minimap_kmer_size}")
    logger.info("-" * 40)
    logger.info("All input parameters validated successfully.")

    pipeline_start_time = time.time()

    # Print summary to console.
    print("-" * 40)
    print("Input Parameters Summary:")
    print("-" * 40)
    params_list = [
        ("Mode", "analyze"),
        ("Input(s)", ', '.join(input_fastq)),
        ("Output directory", out_directory),
        ("Reference", input_reference),
        ("Sequencing Type", seq_class)
    ]
    if seq_summary:
        params_list.append(("Sequencing Summary", seq_summary))
    for name, value in params_list:
        print(f"{name}: {value}")
    print("-" * 40)
    print(f"sequenoscope analyze version {__version__}: Analyzing reads...")
    print("-" * 40)

    logger.info("Creating Sequence object for input FASTQ files.")
    sequencing_sample = Sequence("Test", input_fastq)

    logger.info("Extracting reads with FastqExtractor.")
    print("-" * 40)
    print("Processing FASTQ file(s)...")
    print("-" * 40)
    extractor_run = FastqExtractor(sequencing_sample, out_prefix=f"{out_prefix}_read_list", out_dir=intermediate_dir)

    if seq_class.upper() == SequenceTypes.paired_end:
        logger.info("Extracting and renaming paired-end reads.")
        extractor_run.extract_paired_reads()
        rename_read_ids_run = FastqPairedEndRenamer(sequencing_sample, extractor_run.result_files["read_list_file"],
                                                    out_prefix=f"{out_prefix}_renamed_reads", out_dir=intermediate_dir)
        rename_read_ids_run.rename()
        logger.info("Renaming complete. Updating sequence object.")
        sequencing_sample = Sequence("Test", rename_read_ids_run.result_files["fastq_file_renamed"])
        extractor_run = FastqExtractor(sequencing_sample, out_prefix=f"{out_prefix}_read_list", out_dir=intermediate_dir)
        extractor_run.alt_extract_paired_reads()
    else:
        logger.info("Extracting single-end reads.")
        extractor_run.extract_single_reads()

    logger.info("Filtering reads with FastP.")
    fastp_run_process = FastPRunner(
        sequencing_sample,
        intermediate_dir,
        f"{out_prefix}_fastp_output",
        qualified_quality_phred=quality_threshold,
        min_read_len=min_len,
        max_read_len=max_len,
        trim_front_bp=trim_front,
        trim_tail_bp=trim_tail,
        report_only=False,
        dedup=False,
        threads=threads
    )
    fastp_run_process.run_fastp()
    logger.info("Read filtering complete.")

    print("-" * 40)
    print("Mapping FASTQ based on provided reference FASTA file...")
    print("-" * 40)
    logger.info("Running Minimap2 for read mapping.")
    sequencing_sample_filtered = Sequence("Test", fastp_run_process.result_files["output_files_fastp"])
    minimap_run_process = Minimap2Runner(
        sequencing_sample_filtered,
        intermediate_dir,
        input_reference,
        f"{out_prefix}_mapped_sam",
        threads=threads,
        kmer_size=minimap_kmer_size
    )
    minimap_run_process.run_minimap2()
    logger.info("Minimap2 mapping complete.")

    logger.info("Converting SAM to BAM and extracting FASTQ from mapped reads.")
    sam_to_bam_process = SamBamProcessor(
        minimap_run_process.result_files["sam_output_file"],
        intermediate_dir,
        input_reference,
        f"{out_prefix}_mapped_bam",
        thread=threads
    )
    sam_to_bam_process.run_samtools_bam()

    bam_to_fastq_process = SamBamProcessor(
        sam_to_bam_process.result_files["bam_output"],
        intermediate_dir,
        input_reference,
        f"{out_prefix}_mapped_fastq",
        thread=threads
    )
    bam_to_fastq_process.run_samtools_fastq()

    print("-" * 40)
    print("Calculating distances...")
    print("-" * 40)
    logger.info("Running Mash analysis.")
    mash_run = MashSketcher(intermediate_dir, out_prefix)
    mash_results = mash_run.run_mash_sketch(fastp_run_process.result_files["output_files_fastp"])
    mash_genome_size = mash_results["Genome Size"]
    mash_coverage = mash_results["Coverage"]
    logger.info(f"Mash complete. Genome Size: {mash_genome_size}, Coverage: {mash_coverage}")

    print("-" * 40)
    print("Creating manifest files...")
    print("-" * 40)
    logger.info("Creating manifest files and summaries.")
    # Create manifest files in the final output directory (outside intermediates)
    if seq_summary is not None and GeneralSeqParser.check_seq_summary(seq_summary):
        logger.info("Using sequencing summary to create manifest.")
        manifest_with_sum_run = SeqManifest(
            out_prefix,
            sam_to_bam_process.result_files["bam_output"],
            f"{out_prefix}_manifest",
            out_dir=out_directory,
            min_coverage=min_cov,
            fastp_fastq=fastp_run_process.result_files["output_files_fastp"],
            read_list=extractor_run.result_files["read_list_file"],
            in_seq_summary=seq_summary
        )
        fastp_file = GeneralSeqParser(fastp_run_process.result_files["json"], "json")
        seq_summary_single_end_run = SeqManifestSummary(
            out_prefix,
            manifest_with_sum_run.bam_obj,
            f"{out_prefix}_manifest_summary",
            out_dir=out_directory,
            genome_size=mash_genome_size,
            coverage=mash_coverage,
            fastp_json_file=fastp_file.parsed_file
        )
        seq_summary_single_end_run.generate_summary()
        logger.info("Manifest and summary creation with sequencing summary complete.")
    else:
        logger.info("No valid sequencing summary provided. Creating manifest using default time bounds.")
        manifest_no_sum_run = SeqManifest(
            out_prefix,
            sam_to_bam_process.result_files["bam_output"],
            f"{out_prefix}_manifest",
            out_dir=out_directory,
            min_coverage=min_cov,
            fastp_fastq=fastp_run_process.result_files["output_files_fastp"],
            read_list=extractor_run.result_files["read_list_file"],
            in_fastq=input_fastq,
            start_time=start_time_default,
            end_time=end_time_default
        )
        fastp_file = GeneralSeqParser(fastp_run_process.result_files["json"], "json")
        if seq_class.upper() == SequenceTypes.paired_end:
            logger.info("Generating summary for paired-end reads without sequencing summary.")
            seq_summary_no_sum_run = SeqManifestSummary(
                out_prefix,
                manifest_no_sum_run.bam_obj,
                f"{out_prefix}_manifest_summary",
                out_dir=out_directory,
                genome_size=mash_genome_size,
                coverage=mash_coverage,
                fastp_json_file=fastp_file.parsed_file,
                paired=True
            )
        else:
            logger.info("Generating summary for single-end reads without sequencing summary.")
            seq_summary_no_sum_run = SeqManifestSummary(
                out_prefix,
                manifest_no_sum_run.bam_obj,
                f"{out_prefix}_manifest_summary",
                out_dir=out_directory,
                genome_size=mash_genome_size,
                coverage=mash_coverage,
                fastp_json_file=fastp_file.parsed_file,
                paired=False
            )
        seq_summary_no_sum_run.generate_summary()
        logger.info("Manifest and summary creation without sequencing summary complete.")

    pipeline_end_time = time.time()
    total_runtime_seconds = pipeline_end_time - pipeline_start_time

    print("-" * 40)
    print("All Done!")
    print(f"Total runtime: {format_time(total_runtime_seconds)}")
    print("-" * 40)
    logger.info("Analysis pipeline completed successfully.")
    logger.info(f"Total runtime: {format_time(total_runtime_seconds)}")
    logger.info("All operations are complete.")


if __name__ == '__main__':
    run()
