import os
import pytest
from pathlib import Path
from seq_manifest_plots import DataUtil, SourceFileTaxonCoveredBarChart, BoxPlot, RatioBarChart, SingleRatioBarChart
from decision_bar_chart import IndependentDecisionStackedBarChart, CumulativeDecisionBarChart 

def test_read_data():
    path = '/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt'
    df = DataUtil.read_data(path)
    assert not df.empty, 'DataFrame should not be empty'

def test_source_file_taxon_covered_bar_chart():
    chart = SourceFileTaxonCoveredBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt')
    chart.generate_plot()
    assert os.path.exists("source_file_taxon_covered_bar_chart.html"), "HTML file not generated"

def test_box_plot():
    plot = BoxPlot('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt', 'taxon_%_covered_bases')
    plot.generate_plot()
    assert os.path.exists("box_plot.html"), "HTML file not generated"

def test_ratio_bar_chart():
    chart = RatioBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt')
    chart.generate_plot()
    assert os.path.exists("ratio_bar_chart.html"), "HTML file not generated"

def test_single_ratio_bar_chart():
    chart = SingleRatioBarChart('/home/ameknas/sequenoscope-1/test_stool_samples/0023/Qia_sample_manifest_summary.txt', '/home/ameknas/sequenoscope-1/test_stool_samples/0023/ZYMO_sample_manifest_summary.txt', 'taxon_%_covered_bases')
    chart.generate_plot()
    assert os.path.exists("single_ratio_bar_chart.html"), "HTML file not generated"

#----------------------------------------

@pytest.fixture
def sample_data_1():
    return 'path_to_sample_file_1.csv'


def test_IndependentDecisionStackedBarChart_process_data(sample_data_1):
    chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='hours')
    chart.process_data()
    assert hasattr(chart, 'total_count_2')
    assert hasattr(chart, 'decision_count')

def test_IndependentDecisionStackedBarChart_create_trace(sample_data_1):
    chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='hours')
    chart.create_trace()
    assert hasattr(chart, 'hourly_counts')
    assert hasattr(chart, 'count_values')

def test_IndependentDecisionStackedBarChart_create_chart(sample_data_1):
    chart = IndependentDecisionStackedBarChart(data_path=sample_data_1, time_bin_unit='hours')
    chart.create_chart()
    chart_file = Path("independent_decision_bar_chart.html")
    assert chart_file.is_file()
    chart_file.unlink()  # delete the file after the test

def test_CumulativeDecisionBarChart_create_chart(sample_data_1):
    chart = CumulativeDecisionBarChart(data_path=sample_data_1, time_bin_unit='hours')
    chart.create_chart()
    chart_file = Path("cumulative_decision_bar_chart.html")
    assert chart_file.is_file()
    chart_file.unlink()