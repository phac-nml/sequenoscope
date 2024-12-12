#!/usr/bin/env python
import os
import sys
import time
import logging
import warnings
import pandas as pd
pd.options.mode.chained_assignment = None
import argparse as ap

from sequenoscope.utils.__init__ import format_time
from sequenoscope.plot.seq_manifest_plots import SeqManifestPlotter
from sequenoscope.plot.decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart
from sequenoscope.plot.stats_table import MakeStatsTable
from sequenoscope.plot.violin_plot import ViolinPlotter
from sequenoscope.version import __version__

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scipy.stats.morestats')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in long_scalars')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in double_scalars')

def parse_args():
    parser = ap.ArgumentParser(
        prog="sequenoscope",
        usage="sequenoscope plot --test_dir <test_dir_path> --control_dir <control_dir_path> --output_dir <out_path>\nFor help use: sequenoscope plot -h or sequenoscope plot --help", 
        description="%(prog)s version {}: a flexible tool for processing multiplatform sequencing data: analyze, subset/filter, compare and visualize.".format(__version__),
        formatter_class= ap.RawTextHelpFormatter
    )
    
    parser._optionals.title = "Optional Arguments"

    # Required Paths Group
    paths_group = parser.add_argument_group('Required Paths', 'Specify the necessary directories for the tool.')
    paths_group.add_argument('-T', '--test_dir', type=str, required=True, help="Path to test directory.\n\n")
    paths_group.add_argument('-C', '--control_dir', type=str, required=True, help="Path to control directory.\n\n")
    paths_group.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory designation.\n\n")
    paths_group.add_argument('--force', action='store_true', help='Force overwrite of existing results directory.\n\n')

    # Plotting Options Group
    plotting_group = parser.add_argument_group('Plotting Options', 'Customize the appearance and data for plots.\n\n')
    plotting_group.add_argument('-op', '--output_prefix', type=str, default='sample', help="Output prefix added before plot names. Default is 'sample'.\n\n")
    plotting_group.add_argument('-AS', '--adaptive_sampling', action="store_true", help="Generate decision bar charts for adaptive sampling if utilized during sequencing.\n\n")
    plotting_group.add_argument('-single', '--single_charts', type=bool, default=False, help="Generate charts for data based on selected comparison metric.\n\n")
    plotting_group.add_argument('--comparison_metric', default='taxon_%_covered_bases', 
                                choices=['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases', 'mean_read_length', 'taxon_length', 'taxon_%_covered_bases', 'taxon_mean_read_length'], 
                                type=str, help='Parameter for the box plot and single ratio bar chart. Default: taxon_%_covered_bases.\n\n')
    plotting_group.add_argument('-VP', '--violin_data_percent', default=0.1, type=float, help='Fraction of the data to use for the violin plot. Default=0.1.\n\n')
    plotting_group.add_argument('-bin', '--time_bin_unit', default="minutes", choices=['seconds', 'minutes', '5m', '15m', 'hours'], type=str, help='Time bin used for decision bar charts.\n\n')
    plotting_group.add_argument('-legend', '--taxon_chart_legend', type=bool, default=False, help="Generate a legend for the source file taxon covered bar chart.\n\n")

    return parser.parse_args()

def run():
    args = parse_args()

    test_dir = args.test_dir
    control_dir = args.control_dir
    output_dir = args.output_dir
    output_prefix = args.output_prefix
    adaptive_sampling = args.adaptive_sampling
    single_charts = args.single_charts
    force = args.force
    summary_comp_parameter = args.comparison_metric
    violin_data_precent = args.violin_data_percent
    time_bin_unit = args.time_bin_unit
    taxon_legend = args.taxon_chart_legend

    # Setup output directory and logging
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir, 0o755)
    elif not force:
        print(f"Error directory {output_dir} already exists, if you want to overwrite existing results then specify --force")
        sys.exit()

    log_filename = "plot.log"
    log_filepath = os.path.join(output_dir, log_filename)
    
    logger = logging.getLogger("sequenoscope_plot")
    logger.setLevel(logging.INFO)
    
    fh = logging.FileHandler(log_filepath, mode='w')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Begin logging
    logger.info("Starting 'sequenoscope plot' module.")
    logger.info(f"Version: {__version__}")
    logger.info("Validating input parameters...")

    # Verify required files are found
    params_list = [
        ("Mode", "plot"),
        ("Test directory", test_dir),
        ("Control directory", control_dir),
        ("Output directory", output_dir),
        ("Adaptive Sampling", adaptive_sampling)
    ]

    logger.info("Input Parameters:")
    logger.info("-" * 40)
    for param_name, param_value in params_list:
        logger.info(f"{param_name}: {param_value}")
    logger.info(f"Comparison metric: {summary_comp_parameter}")
    logger.info(f"Violin data fraction: {violin_data_precent}")
    logger.info(f"Time bin unit: {time_bin_unit}")
    logger.info(f"Taxon chart legend: {taxon_legend}")
    logger.info(f"Single charts: {single_charts}")
    logger.info("-" * 40)
    logger.info("All input parameters validated successfully.")

    print("-" * 40)
    print("Input Parameters Summary:")
    print("-" * 40)
    for param_name, param_value in params_list:
        print(f"{param_name}: {param_value}")
    print("-"*40)
    print(f"sequenoscope plot version {__version__}: Extracting files...")
    print("-"*40)

    start_time = time.time()

    logger.info("Searching for required manifest and summary files in specified directories.")
    test_manifest = None
    test_manifest_summary = None
    control_manifest = None
    control_manifest_summary = None

    for f in os.listdir(test_dir):
        if 'manifest.txt' in f:
            test_manifest = os.path.join(test_dir, f)
        elif 'manifest_summary.txt' in f:
            test_manifest_summary = os.path.join(test_dir, f)

    for f in os.listdir(control_dir):
        if 'manifest.txt' in f:
            control_manifest = os.path.join(control_dir, f)
        elif 'manifest_summary.txt' in f:
            control_manifest_summary = os.path.join(control_dir, f)

    if not test_manifest or not control_manifest or not test_manifest_summary or not control_manifest_summary:
        logger.error("One or more required files (manifest or manifest_summary) were not found.")
        raise ValueError("One or more required files were not found in the given directories.")

    logger.info("Required files found successfully.")

    print("-"*40)
    print("Plotting manifest summary plots...")
    print("-"*40)
    logger.info("Generating taxon covered bar charts, ratio charts, and stats tables.")

    taxon_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    ratio_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    stats_table = MakeStatsTable(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)

    taxon_bar_chart.generate_source_file_taxon_covered_bar_chart(show_legend=taxon_legend)
    ratio_bar_chart.generate_ratio_bar_chart()

    stats_table.generate_stats()
    stats_table.save_to_csv()
    logger.info("Taxon bar charts, ratio charts, and stats tables generated successfully.")

    if single_charts:
        logger.info("Generating single charts based on the specified comparison metric.")
        box_plot = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
        single_ratio_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)

        box_plot.generate_box_plot(summary_comp_parameter)
        single_ratio_bar_chart.generate_single_ratio_bar_chart(summary_comp_parameter)
        logger.info("Single box plot and single ratio bar charts generated successfully.")

    print("-"*40)
    print("Plotting seq manifest plots...")
    print("-"*40)
    logger.info("Generating violin plots for read quality score and read length.")

    violin_plot_read_qscore = ViolinPlotter(test_manifest, control_manifest, output_dir, output_prefix,
                                            quality_metric='read_qscore', fraction=violin_data_precent)
    violin_plot_read_length = ViolinPlotter(test_manifest, control_manifest, output_dir, output_prefix,
                                            quality_metric='read_len', fraction=violin_data_precent)
        
    violin_plot_read_qscore.generate_chart()
    violin_plot_read_length.generate_chart()
    logger.info("Violin plots generated successfully.")

    if adaptive_sampling:
        logger.info("Adaptive sampling detected. Generating decision bar charts.")
        test_independent_decision_bar_chart = IndependentDecisionStackedBarChart(test_manifest, output_dir, output_prefix + "_test", time_bin_unit)
        control_independent_decision_bar_chart = IndependentDecisionStackedBarChart(control_manifest, output_dir, output_prefix + "_control", time_bin_unit)

        test_cumulative_decision_bar_chart = CumulativeDecisionBarChart(test_manifest, output_dir, output_prefix + "_test", time_bin_unit)
        control_cumulative_decision_bar_chart = CumulativeDecisionBarChart(control_manifest, output_dir, output_prefix + "_control", time_bin_unit)

        test_independent_decision_bar_chart.generate_chart()
        control_independent_decision_bar_chart.generate_chart()
        test_cumulative_decision_bar_chart.generate_chart()
        control_cumulative_decision_bar_chart.generate_chart()
        logger.info("Decision bar charts generated successfully.")

    end_time = time.time()
    total_runtime_minutes = end_time - start_time

    logger.info("Plotting pipeline completed successfully.")
    logger.info(f"Total runtime: {format_time(total_runtime_minutes)}")
    logger.info("All operations are complete.")

    print("-"*40)
    print("All Done!")
    print(f"total runtime: {format_time(total_runtime_minutes)}")
    print("-"*40)
