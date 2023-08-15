#!/usr/bin/env python
import os
import sys
import argparse as ap
from sequenoscope.plot.seq_manifest_plots import SeqManifestPlotter
from sequenoscope.plot.decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart
from sequenoscope.plot.stats_table import MakeStatsTable
from sequenoscope.plot.violin_plot import ViolinPlotter
from sequenoscope.version import __version__

# TODO: add comments to all classes to describe what's going on
# TODO: Make sure that plot.py has is complete in terms of update and help messages
# TODO: change the errors that pop up when plotting decision bar graphs 



def parse_args():
    parser = ap.ArgumentParser(prog="sequenoscope",
                               usage="sequenoscope plot --test_dir <test_dir_path> --control_dir <control_dir_path> --output_dir <out_path>\nFor help use: sequenoscope plot -h or sequenoscope plot --help", 
                                description="%(prog)s version {}: a tool for analyzing and processing sequencing data.".format(__version__),
                                formatter_class= ap.RawTextHelpFormatter)
    
    # Customize the title of the optional arguments section in the help menu
    parser._optionals.title = "Optional Arguments"

    # Required Paths Group
    paths_group = parser.add_argument_group('Required Paths', 'Specify the necessary directories for the tool.')
    paths_group.add_argument('-T', '--test_dir', type=str, required=True, help="Path to test directory.\n\n")
    paths_group.add_argument('-C', '--control_dir', type=str, required=True, help="Path to control directory.\n\n")
    paths_group.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory designation.\n\n")
    paths_group.add_argument('--force', action='store_true', help='Force overwrite of existing results directory.\n\n')
    # Plotting Options Group
    plotting_group = parser.add_argument_group('Plotting Options', 'Customize the appearance and data for plots.\n\n')
    plotting_group.add_argument('-o_pre', '--output_prefix', type=str, default='sample', help="Output prefix added before plot names. Default is 'sample'.\n\n")
    plotting_group.add_argument('--comparison_metric', default='taxon_%_covered_bases', choices=['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases', 'mean_read_length', 'taxon_length', 'taxon_covered_bases', 'taxon_%_covered_bases', 'taxon_mean_read_length'], type=str, help='Type of parameter for the box plot and single ratio bar chart. Default parameter is taxon_%%_covered_bases.\n\n')
    plotting_group.add_argument('-VP', '--violin_data_percent', default=0.1, type=float, help='Fraction of the data to use for the violin plot.\n\n')
    plotting_group.add_argument('-bin', '--time_bin_unit', default="minutes", choices=['seconds', 'minutes', '5m', '15m', 'hours'], type=str, help='Time bin used for decision bar charts.\n\n')

    return parser.parse_args()

def run():
    args = parse_args()

    test_dir = args.test_dir
    control_dir = args.control_dir
    output_dir = args.output_dir
    output_prefix = args.output_prefix
    force = args.force
    summary_comp_parameter = args.comparison_metric
    violin_data_precent = args.violin_data_percent
    time_bin_unit = args.time_bin_unit

    print("-"*40)
    print(f"sequenoscope plot version {__version__}: Extracting files from directories for plotting")
    print("-"*40)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir, 0o755)
    elif not force:
        print(f"Error directory {output_dir} already exists, if you want to overwrite existing results then specify --force")
        sys.exit()

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
        raise ValueError("One or more required files were not found in the given directories.")
    
    print("-"*40)
    print("Plotting manifest summary plots")
    print("-"*40)

    taxon_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    box_plot = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    ratio_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    single_ratio_bar_chart = SeqManifestPlotter(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)
    stats_table = MakeStatsTable(test_manifest_summary, control_manifest_summary, output_dir, output_prefix=output_prefix)

    taxon_bar_chart.generate_source_file_taxon_covered_bar_chart()
    box_plot.generate_box_plot(summary_comp_parameter)
    ratio_bar_chart.generate_ratio_bar_chart()
    single_ratio_bar_chart.generate_single_ratio_bar_chart(summary_comp_parameter)
    stats_table.generate_stats()
    stats_table.save_to_csv()

    print("-"*40)
    print("Plotting seq manifest plots")
    print("-"*40)

    test_independent_decision_bar_chart = IndependentDecisionStackedBarChart(test_manifest, output_dir, "test_" + output_prefix, time_bin_unit)
    control_independent_decision_bar_chart = IndependentDecisionStackedBarChart(control_manifest, output_dir, "control_" + output_prefix, time_bin_unit)
    test_cumulative_decision_bar_chart = CumulativeDecisionBarChart(test_manifest, output_dir, "test_" + output_prefix, time_bin_unit)
    control_cumulative_decision_bar_chart = CumulativeDecisionBarChart(control_manifest, output_dir, "control_" + output_prefix, time_bin_unit)
    violin_plot_read_qscore = ViolinPlotter(test_manifest, control_manifest, output_dir, "read_qscore_" + output_prefix, quality_metric='read_qscore', fraction=violin_data_precent)
    violin_plot_read_length = ViolinPlotter(test_manifest, control_manifest, output_dir, "read_len_" + output_prefix, quality_metric='read_len', fraction=violin_data_precent)


    test_independent_decision_bar_chart.generate_chart()
    control_independent_decision_bar_chart.generate_chart()
    test_cumulative_decision_bar_chart.generate_chart()
    control_cumulative_decision_bar_chart.generate_chart()
    violin_plot_read_qscore.generate_chart()
    violin_plot_read_length.generate_chart()

    print("-"*40)
    print("All Done")
    print("-"*40)