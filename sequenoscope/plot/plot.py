#!/usr/bin/env python
import os
import sys
import time
import logging
import warnings
import pandas as pd
import argparse as ap

from sequenoscope.utils.__init__ import format_time
from sequenoscope.plot.seq_manifest_plots import SeqManifestPlotter
from sequenoscope.plot.summary_table import SummaryTable
from sequenoscope.plot.violin_plot import ViolinPlotter
from sequenoscope.plot.decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart
from sequenoscope.version import __version__

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scipy.stats.morestats')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in long_scalars')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in double_scalars')

def parse_args():
    parser = ap.ArgumentParser(
        prog="sequenoscope",
        usage="sequenoscope plot --test_dir <test_dir_path> --control_dir <control_dir_path> --output_dir <out_path>\nFor help use: sequenoscope plot -h or sequenoscope plot --help",
        description=f"sequenoscope version {__version__}: a flexible tool for processing multiplatform sequencing data: analyze, subset/filter, compare and visualize.",
        formatter_class=ap.RawTextHelpFormatter
    )
    
    parser._optionals.title = "Optional Arguments"

    # Required Paths Group
    paths_group = parser.add_argument_group('Required Paths', 'Specify the necessary directories for the tool.')
    paths_group.add_argument('-T', '--test_dir', type=str, required=True, help="Path to test directory.\n")
    paths_group.add_argument('-C', '--control_dir', type=str, required=True, help="Path to control directory.\n")
    paths_group.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory designation.\n")
    paths_group.add_argument('--force', action='store_true', help='Force overwrite of existing results directory.\n')

    # Plotting Options Group
    plotting_group = parser.add_argument_group('Plotting Options', 'Customize the appearance and data for plots.\n')
    plotting_group.add_argument('-op', '--output_prefix', type=str, default='sample', help="Output prefix added before plot names. Default is 'sample'.\n")
    plotting_group.add_argument('-AS', '--adaptive_sampling', action="store_true", help="Generate decision bar charts for adaptive sampling if utilized during sequencing.\n")
    plotting_group.add_argument('-VP', '--violin_data_percent', default=0.1, type=float, help='Fraction of the data to use for the violin plot. Default=0.1.\n')
    plotting_group.add_argument('-bin', '--time_bin_unit', default="minutes", choices=['seconds', 'minutes', '5m', '15m', 'hours'], type=str, help='Time bin used for decision bar charts.\n')

    return parser.parse_args()

def run():
    args = parse_args()

    test_dir = args.test_dir
    control_dir = args.control_dir
    output_dir = args.output_dir
    output_prefix = args.output_prefix
    adaptive_sampling = args.adaptive_sampling
    force = args.force
    violin_data_percent = args.violin_data_percent
    time_bin_unit = args.time_bin_unit

    # Setup output directory and logging
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir, 0o755)
    elif not force:
        print(f"Error: Directory {output_dir} already exists. Use --force to overwrite.")
        sys.exit()

    log_filepath = os.path.join(output_dir, "plot.log")
    logger = logging.getLogger("sequenoscope_plot")
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_filepath, mode='w')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("Starting 'sequenoscope plot' module.")
    logger.info(f"Version: {__version__}")
    logger.info("Validating input parameters...")

    params_list = [
        ("Mode", "plot"),
        ("Test directory", test_dir),
        ("Control directory", control_dir),
        ("Output directory", output_dir),
        ("Adaptive Sampling", adaptive_sampling)
    ]
    logger.info("Input Parameters:")
    logger.info("-" * 40)
    for name, value in params_list:
        logger.info(f"{name}: {value}")
    logger.info(f"Violin data fraction: {violin_data_percent}")
    logger.info(f"Time bin unit: {time_bin_unit}")
    logger.info("-" * 40)

    print("-" * 40)
    print("Input Parameters Summary:")
    for name, value in params_list:
        print(f"{name}: {value}")
    print("-" * 40)
    print(f"sequenoscope plot version {__version__}: Extracting files...")
    print("-" * 40)

    start_time = time.time()
    logger.info("Searching for required manifest and summary files in specified directories.")

    test_manifest = None
    test_manifest_summary = None
    control_manifest = None
    control_manifest_summary = None

    for f in os.listdir(test_dir):
        if 'manifest.txt' in f and 'summary' not in f:
            test_manifest = os.path.join(test_dir, f)
        elif 'manifest_summary.txt' in f:
            test_manifest_summary = os.path.join(test_dir, f)

    for f in os.listdir(control_dir):
        if 'manifest.txt' in f and 'summary' not in f:
            control_manifest = os.path.join(control_dir, f)
        elif 'manifest_summary.txt' in f:
            control_manifest_summary = os.path.join(control_dir, f)

    if not all([test_manifest, control_manifest, test_manifest_summary, control_manifest_summary]):
        logger.error("One or more required files (manifest or manifest_summary) were not found.")
        raise ValueError("Required files not found in the provided directories.")

    logger.info("Required files found successfully.")

    print("-" * 40)
    print("Plotting manifest summary plots...")
    print("-" * 40)
    logger.info("Generating taxon covered bar charts and summary table.")

    # Generate taxon covered bar chart using summary files, always show legend.
    taxon_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    taxon_bar_chart.generate_source_file_taxon_covered_bar_chart()

    # Generate summary table
    from sequenoscope.plot.summary_table import SummaryTable  # ensure the correct module is imported
    summary_table = SummaryTable(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    summary_table.generate_summary()
    summary_table.save_to_csv()
    logger.info("Taxon bar charts and summary table generated successfully.")

    print("-" * 40)
    print("Plotting seq manifest plots...")
    print("-" * 40)
    logger.info("Generating violin plots for read quality score and read length.")

    violin_plot_read_qscore = ViolinPlotter(test_manifest, control_manifest, output_dir, output_prefix,
                                             quality_metric='read_qscore', fraction=violin_data_percent)
    violin_plot_read_length = ViolinPlotter(test_manifest, control_manifest, output_dir, output_prefix,
                                            quality_metric='read_len', fraction=violin_data_percent)
        
    violin_plot_read_qscore.generate_chart()
    violin_plot_read_length.generate_chart()
    logger.info("Violin plots generated successfully.")

    if adaptive_sampling:
        logger.info("Adaptive sampling detected. Generating decision bar charts.")
        # Create a subdirectory for decision bar charts that includes the time bin unit in its name
        decision_bar_dir = os.path.join(output_dir, f"decision_bar_charts_{time_bin_unit}")
        if not os.path.exists(decision_bar_dir):
            os.mkdir(decision_bar_dir, 0o755)
        test_independent_chart = IndependentDecisionStackedBarChart(test_manifest, decision_bar_dir, output_prefix + "_test", time_bin_unit)
        control_independent_chart = IndependentDecisionStackedBarChart(control_manifest, decision_bar_dir, output_prefix + "_control", time_bin_unit)
        test_cumulative_chart = CumulativeDecisionBarChart(test_manifest, decision_bar_dir, output_prefix + "_test", time_bin_unit)
        control_cumulative_chart = CumulativeDecisionBarChart(control_manifest, decision_bar_dir, output_prefix + "_control", time_bin_unit)
        test_independent_chart.generate_chart()
        control_independent_chart.generate_chart()
        test_cumulative_chart.generate_chart()
        control_cumulative_chart.generate_chart()
        logger.info("Decision bar charts generated successfully.")

    # Generate default chart comparing taxon mean read length for each species
    logger.info("Generating default chart comparing taxon mean read length.")
    default_mean_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    default_mean_chart.generate_mean_read_length_chart()

    # Generate default chart comparing taxon mean coverage for each species
    logger.info("Generating default chart comparing taxon mean coverage.")
    default_coverage_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    default_coverage_chart.generate_mean_coverage_chart()

    end_time = time.time()
    total_runtime = end_time - start_time
    logger.info("Plotting pipeline completed successfully.")
    logger.info(f"Total runtime: {format_time(total_runtime)}")

    print("-" * 40)
    print("All Done!")
    print(f"Total runtime: {format_time(total_runtime)}")
    print("-" * 40)

if __name__ == '__main__':
    run()
