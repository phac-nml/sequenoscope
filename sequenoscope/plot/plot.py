#!/usr/bin/env python
import os
import sys
import argparse as ap
from sequenoscope.plot.seq_manifest_plots import SeqManifestPlotter
from sequenoscope.plot.decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart
from sequenoscope.plot.stats_table import MakeStatsTable
from sequenoscope.plot.violin_plot import ViolinPlotter
from sequenoscope.version import __version__

def parse_args():
    parser = ap.ArgumentParser(prog="sequenoscope",
                               usage="sequenoscope plot <'manifest', 'summary'> --test_file <test_file_path> --control_file <control_file_path> --plot_type <plot_type> --output_file <out_path>\nFor help use: sequenoscope plot -h or sequenoscope plot --help", 
                                description="%(prog)s version {}: a tool for analyzing and processing sequencing data.".format(__version__))
                                # formatter_class= ap.RawTextHelpFormatter)


    parser._optionals.title = "Arguments"

    subparsers = parser.add_subparsers(dest='plot_type')
    subparsers.required = True

    parser_common = ap.ArgumentParser(add_help=False)
    parser_common.add_argument('--output_file', type=str, required=True, help="[REQUIRED] Output directory designation")
    parser_common.add_argument('--force', required=False, help='Force overwite of existing results directory', action='store_true')

    # Sub-parser for plot types that require manifest files
    summary_parser = subparsers.add_parser('summary', parents=[parser_common], help='yp')
    summary_parser._optionals.title = "Arguments"
    summary_parser.add_argument('--test_file', type=str, required=True, help="[REQUIRED] Path to test file to process.")
    summary_parser.add_argument('--control_file', type=str, required=True, help="[REQUIRED] Path to control file to process.")
    summary_parser.add_argument('--plot_type', type=str, default= 'all', choices=['taxon_bar_chart', 'boxplot', 'ratio_bar_chart', 'single_ratio_bar_chart', 'stats_table', 'all'], required=True, help='the type of plot to generate. default is that all plots are generated')
    summary_parser.add_argument('--comp_parameter', default = 'taxon_%%_covered_bases', choices=['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases', 'mean_read_length', 'taxon_length',	'taxon_covered_bases', 'taxon_%%_covered_bases', 'taxon_mean_read_length'], type=str, help='type of paramter for the box plot and single ratio bar chart. default paramter is taxon_%%_covered_bases')
    
    # Sub-parser for plot types that require summary files
    manifest_parser = subparsers.add_parser('manifest', parents=[parser_common])
    manifest_parser.add_argument('--test_file', type=str, required=True, help="[REQUIRED] Path to test file to process.")
    manifest_parser.add_argument('--time_bin_unit', default="minutes", choices=['seconds', 'minutes', 'hours'], type=str, help='time bin used for decision bar charts')
    manifest_parser.add_argument('--control_file', type=str, required=True, help="[REQUIRED] Path to control file to process.")
    manifest_parser.add_argument('--plot_type', default='all', choices=['IND_decision_bar_chart', 'CUM_decision_bar_chart', 'violin_plot', 'all'], help='the type of plot to generate. default is that all plots are generated')
    manifest_parser.add_argument('--fraction', default = 0.1, type=float, help='fraction of the data to use for the violin plot')
    manifest_parser.add_argument('--violin_parameter', default = 'read_qscore', choices=['read_qscore', 'read_len'], type=str, help='paramter used for generating the violin plot. default is read_qscore.')

    return parser.parse_args()

def run():
    args = parse_args()

    plot_type = args.plot_type
    output_file = args.output_file
    force = args.force

    if plot_type == 'summary':
        test_file_summary = args.test_file
        control_file_summary = args.control_file
        plot_type_summary = args.plot_type
        comp_parameter = args.comp_parameter

    elif plot_type == 'manifest':
        test_file_manifest = args.test_file
        control_file_manifest = args.control_file
        time_bin_unit = args.time_bin_unit
        plot_type_manifest = args.plot_type
        fraction = args.fraction
        violin_parameter = args.violin_parameter

    print("-"*40)
    print(f"sequenoscope plot version {__version__}: producing the plots based on input files")
    print("-"*40)

    ## intializing directory for files

    if not os.path.isdir(output_file):
        os.mkdir(output_file, 0o755)
    elif not force:
        print(f"Error directory {output_file} already exists, if you want to overwrite existing results then specify --force")
        sys.exit()

    if plot_type == 'summary':
        taxon_bar_chart = SeqManifestPlotter(test_file_summary, control_file_summary)
        box_plot = SeqManifestPlotter(test_file_summary, control_file_summary)
        ratio_bar_chart = SeqManifestPlotter(test_file_summary, control_file_summary)
        single_ratio_bar_chart = SeqManifestPlotter(test_file_summary, control_file_summary)
        stats_table = MakeStatsTable(test_file_summary, control_file_summary)
        if plot_type_summary == 'all':
            taxon_bar_chart.generate_source_file_taxon_covered_bar_chart()
            box_plot.generate_box_plot(comp_parameter)
            ratio_bar_chart.generate_ratio_bar_chart()
            single_ratio_bar_chart.generate_single_ratio_bar_chart(comp_parameter)
            stats_table.generate_stats()
            stats_table.save_to_csv()
        elif plot_type_summary == 'taxon_bar_chart':
            taxon_bar_chart.generate_source_file_taxon_covered_bar_chart()
        elif plot_type_summary == 'boxplot':
            box_plot.generate_box_plot(comp_parameter)
        elif plot_type_summary == 'ratio_bar_chart':
            ratio_bar_chart.generate_ratio_bar_chart()
        elif plot_type_summary == 'single_ratio_bar_chart':
            single_ratio_bar_chart.generate_single_ratio_bar_chart(comp_parameter)
        elif plot_type_summary == 'stats_table':
            stats_table.generate_stats()
            stats_table.save_to_csv()

    if plot_type == 'manifest':
        test_independent_decision_bar_chart = IndependentDecisionStackedBarChart(test_file_manifest, time_bin_unit)
        control_independent_decision_bar_chart = IndependentDecisionStackedBarChart(control_file_manifest, time_bin_unit)
        test_cumulative_decision_bar_chart = CumulativeDecisionBarChart(test_file_manifest, time_bin_unit)
        control_cumulative_decision_bar_chart = CumulativeDecisionBarChart(control_file_manifest, time_bin_unit)
        violin_plot = ViolinPlotter(test_file_manifest, control_file_manifest, violin_parameter, fraction)
        if plot_type_manifest == 'all':
            test_independent_decision_bar_chart.process_data()
            test_independent_decision_bar_chart.create_trace()
            test_independent_decision_bar_chart.create_chart()
            #test_independent_decision_bar_chart.song_and_dance()

            control_independent_decision_bar_chart.process_data()
            control_independent_decision_bar_chart.create_trace()
            control_independent_decision_bar_chart.create_chart()

            test_cumulative_decision_bar_chart.process_data()
            test_cumulative_decision_bar_chart.create_trace()
            test_cumulative_decision_bar_chart.create_chart()

            control_cumulative_decision_bar_chart.process_data()
            control_cumulative_decision_bar_chart.create_trace()
            control_cumulative_decision_bar_chart.create_chart()

            violin_plot.process_files()
            violin_plot.create_violin_plot()
        elif plot_type_manifest == 'IND_decision_bar_chart':
            test_independent_decision_bar_chart.process_data()
            test_independent_decision_bar_chart.create_trace()
            test_independent_decision_bar_chart.create_chart()

            control_independent_decision_bar_chart.process_data()
            control_independent_decision_bar_chart.create_trace()
            control_independent_decision_bar_chart.create_chart()
        elif plot_type_manifest == 'CUM_decision_bar_chart':
            test_cumulative_decision_bar_chart.process_data()
            test_cumulative_decision_bar_chart.create_trace()
            test_cumulative_decision_bar_chart.create_chart()

            control_cumulative_decision_bar_chart.process_data()
            control_cumulative_decision_bar_chart.create_trace()
            control_cumulative_decision_bar_chart.create_chart()
        elif plot_type_manifest == 'violin_plot':
            violin_plot.process_files()
            violin_plot.create_violin_plot()
    pass

