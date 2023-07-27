import os
import pytest
from pathlib import Path
from seq_manifest_plots import DataUtil, SourceFileTaxonCoveredBarChart, BoxPlot, RatioBarChart, SingleRatioBarChart
from decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart 
from violin_plot import ViolinPlotter
from stats_table import MakeStatsTable

# def test_read_data():
#     path = '/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt'
#     df = DataUtil.read_data(path)
#     assert not df.empty, 'DataFrame should not be empty'

# def test_source_file_taxon_covered_bar_chart():
#     chart = SourceFileTaxonCoveredBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt')
#     chart.generate_plot()
#     assert os.path.exists("source_file_taxon_covered_bar_chart.html"), "HTML file not generated"

# def test_box_plot():
#     plot = BoxPlot('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt', 'taxon_%_covered_bases')
#     plot.generate_plot()
#     assert os.path.exists("box_plot.html"), "HTML file not generated"

# def test_ratio_bar_chart():
#     chart = RatioBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt')
#     chart.generate_plot()
#     assert os.path.exists("ratio_bar_chart.html"), "HTML file not generated"

# def test_single_ratio_bar_chart():
#     chart = SingleRatioBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt', 'taxon_%_covered_bases')
#     chart.generate_plot()
#     assert os.path.exists("single_ratio_bar_chart.html"), "HTML file not generated"

#----------------------------------------

# @pytest.fixture
# def sample_data_1():
#     return '/home/ameknas/sequenoscope-1/test_SE/sample_manifest.txt'


# def test_IndependentDecisionStackedBarChart_process_data(sample_data_1):
#     chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='minutes')
#     chart.process_data()
#     assert hasattr(chart, 'total_count_2')
#     assert hasattr(chart, 'decision_count')

# def test_IndependentDecisionStackedBarChart_create_trace(sample_data_1):
#     chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='minutes')
#     chart.create_trace()
#     assert hasattr(chart, 'hourly_counts')
#     assert hasattr(chart, 'count_values')

# def test_IndependentDecisionStackedBarChart_create_chart(sample_data_1):
#     chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='minutes')
#     chart.create_chart()
#     chart_file = Path("independent_decision_bar_chart.html")
#     assert chart_file.is_file()
#     chart_file.unlink()  # delete the file after the test

# def test_CumulativeDecisionBarChart_create_chart(sample_data_1):
#     chart = CumulativeDecisionBarChart(data_path=sample_data_1, time_bin_unit='minutes')
#     chart.create_chart()
#     chart_file = Path("cumulative_decision_bar_chart.html")
#     assert chart_file.is_file()
#     chart_file.unlink()

#----------------------------------------------------------------------

# def test_violin_plotter():
#     test_file = '/home/ameknas/sequenoscope-1/test_SE/sample_manifest.txt'
#     control_file = '/home/ameknas/sequenoscope-1/test_PE/sample_manifest.txt'
#     plotter = ViolinPlotter(test_file, control_file)

#     # Test initial attributes
#     assert plotter.test_file == test_file
#     assert plotter.control_file == control_file
#     assert plotter.quality_metric == 'read_qscore'
#     assert plotter.fraction == 0.1
#     assert plotter.data is None

#     # Test process_files method
#     plotter.process_files()
#     assert plotter.data is not None

#     # Test create_violin_plot method
#     plotter.create_violin_plot()
#     assert os.path.isfile('violin_comparison_plot.html')

#     # Test if ValueError is raised when quality_metric column doesn't exist
#     with pytest.raises(ValueError):
#         invalid_plotter = ViolinPlotter(test_file, control_file, 'invalid_column')
#         invalid_plotter.process_files()

#     # Test if ValueError is raised when trying to plot before processing files
#     with pytest.raises(ValueError):
#         invalid_plotter = ViolinPlotter(test_file, control_file)
#         invalid_plotter.create_violin_plot()

def test_MakeStatsTable():
    # Define the paths to your test files
    test_file_path = "/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt"
    control_file_path = "/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt"

    stats_table = MakeStatsTable(test_file_path, control_file_path)
    stats_table.generate_stats()
    stats_table.save_to_csv()

    # Check if the output file was created
    assert os.path.isfile('result.csv')