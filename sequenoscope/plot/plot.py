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
                                description="%(prog)s version {}: a tool for analyzing and processing sequencing data.".format(__version__),
                                formatter_class= ap.RawTextHelpFormatter)
    
    parser._optionals.title = "Arguments"
    
    # Main arguments
    parser.add_argument('-T', '--test_dir', type=str, required=True, help="[REQUIRED] Path to test directory.")
    parser.add_argument('-C', '--control_dir', type=str, required=True, help="[REQUIRED] Path to control directory.")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="[REQUIRED] Output directory designation.")
    parser.add_argument('-o_pre', '--output_prefix', type=str, default = 'sample', help="[REQUIRED] Output prefix added before plot names. default is sample")
    parser.add_argument('--force', required=False, help='Force overwrite of existing results directory', action='store_true')
    parser.add_argument('-SCP', '--summary_comp_paramter', default = 'taxon_%_covered_bases', choices=['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases', 'mean_read_length', 'taxon_length',	'taxon_covered_bases', 'taxon_%_covered_bases', 'taxon_mean_read_length'], type=str, help='type of paramter for the box plot and single ratio bar chart. default paramter is taxon_%_covered_bases')
    parser.add_argument('-VP', '--violin_data_percent', default = 0.1, type=float, help='fraction of the data to use for the violin plot')
    parser.add_argument('-bin', '--time_bin_unit', default="minutes", choices=['seconds', 'minutes', 'hours'], type=str, help='time bin used for decision bar charts')

    return parser.parse_args()

def run():
    args = parse_args()

    test_dir = args.test_dir
    control_dir = args.control_dir
    output_dir = args.output_dir
    output_prefix = args.output_prefix
    force = args.force
    summary_comp_parameter = args.summary_comp_paramter
    violin_data_precent = args.violin_data_percent
    time_bin_unit = args.time_bin_unit

    print("-"*40)
    print(f"sequenoscope plot version {__version__}: Extracting files from directories for plotting")
    print("-"*40)


