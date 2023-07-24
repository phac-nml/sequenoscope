import os
import pytest
from seq_manifest_plots import DataUtil, SourceFileTaxonCoveredBarChart, BoxPlot, RatioBarChart, SingleRatioBarChart

def test_read_data():
    path = 'test_file.csv'
    df = DataUtil.read_data(path)
    assert not df.empty, 'DataFrame should not be empty'

def test_source_file_taxon_covered_bar_chart():
    chart = SourceFileTaxonCoveredBarChart('test_file.csv', 'control_file.csv')
    chart.generate_plot()
    assert os.path.exists("source_file_taxon_covered_bar_chart.html"), "HTML file not generated"

def test_box_plot():
    plot = BoxPlot('test_file.csv', 'control_file.csv', 'some_param')
    plot.generate_plot()
    assert os.path.exists("box_plot.html"), "HTML file not generated"

def test_ratio_bar_chart():
    chart = RatioBarChart('test_file.csv', 'control_file.csv')
    chart.generate_plot()
    assert os.path.exists("ratio_bar_chart.html"), "HTML file not generated"

def test_single_ratio_bar_chart():
    chart = SingleRatioBarChart('test_file.csv', 'control_file.csv', 'some_param')
    chart.generate_plot()
    assert os.path.exists("single_ratio_bar_chart.html"), "HTML file not generated"
