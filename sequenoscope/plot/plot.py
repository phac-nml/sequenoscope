#!/usr/bin/env python
import argparse as ap
from sequenoscope.plot.seq_manifest_plots import SeqManifestPlotter
from sequenoscope.plot.decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart
from sequenoscope.plot.stats_table import MakeStatsTable
from sequenoscope.plot.violin_plot import ViolinPlotter
from sequenoscope.version import __version__

def parse_args():
    parser = ap.ArgumentParser(prog="sequenoscope",
                               usage="sequenoscope plot --file_type <'manifest', 'summary'> --test_file <test_file_path> --control_file <control_file_path> --plot-type <plot_type> --output_file <out_path>\nFor help use: sequenoscope plot -h or sequenoscope plot --help", 
                                description="%(prog)s version {}: a tool for analyzing and processing sequencing data.".format(__version__), 
                                formatter_class= ap.RawTextHelpFormatter)

    parser._optionals.title = "Arguments"

    subparsers = parser.add_subparsers(dest='plot_type')
    subparsers.required = True

    parser_common = ap.ArgumentParser(add_help=False)
    parser_common.add_argument('--output-file', type=str, required=True)

    # Sub-parser for plot types that require manifest files
    manifest_parser = subparsers.add_parser('summary', parents=[parser_common])
    manifest_parser.add_argument('--test-file', type=str, required=True)
    manifest_parser.add_argument('--control-file', type=str, required=True)
    manifest_parser.add_argument('--plot-type', default= 'all', choices=['taxon_bar_chart', 'boxplot', 'ratio_bar_chart', 'single_ratio_bar_chart', 'stats_table', 'all'], required=True)
    manifest_parser.add_argument('--comp_parameter', default = 'taxon_%_covered_bases', choices=['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases', 'mean_read_length', 'taxon_length',	'taxon_covered_bases', 'taxon_%_covered_bases', 'taxon_mean_read_length'], type=str)
    
    # Sub-parser for plot types that require summary files
    summary_parser = subparsers.add_parser('manifest', parents=[parser_common])
    summary_parser.add_argument('--test-file', type=str, required=True)
    summary_parser.add_argument('--control-file', type=str, required=True)
    summary_parser.add_argument('--plot-type', default='all', choices=['IND_decision_bar_chart', 'CUM_decision_bar_chart', 'violin_plot', 'all'], required=True)
    summary_parser.add_argument('--fraction', default = 0.1, type=str)
    summary_parser.add_argument('--violin_parameter', default = 'read_qscore', choices=['read_qscore', 'read_len'], type=str)

    return parser.parse_args()

def run():
    args = parse_args()

    plot_type = args.plot_type
    output_file = args.output_file

    if plot_type == 'summary':
        test_file_summary = args.test_file
        control_file_summary = args.control_file
        plot_type_summary = args.plot_type
        comp_parameter = args.comp_parameter

    elif plot_type == 'manifest':
        test_file_manifest = args.test_file
        control_file_manifest = args.control_file
        plot_type_manifest = args.plot_type
        fraction = args.fraction
        violin_parameter = args.violin_parameter
    
    pass

