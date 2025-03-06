#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go

class SeqManifestPlotter:
    def __init__(self, test_file_path, control_file_path, output_dir, output_prefix="sample"):
        """
        Initialize with file paths and output details.
        """
        self.test_file_path = test_file_path
        self.control_file_path = control_file_path
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.color_scale = [
            '#FF7F0E', '#1F77B4', '#FFC0CB', '#2CA02C', '#D62728',
            '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22',
            '#17BECF', '#1A55FF', '#FF991A', '#B2912F', '#197319',
            '#DB2727', '#9C27B0', '#03A9F4', '#FFEB3B', '#009688'
        ]

    def read_data_csv(self, path):
        """
        Read a CSV file and return a DataFrame.
        """
        try:
            return pd.read_csv(path, delimiter='\t')
        except Exception as e:
            print(f"Error reading {path}: {e}")
            return pd.DataFrame()

    def read_and_append_source(self, path, source_name):
        """
        Read CSV data and append a source column.
        """
        df = self.read_data_csv(path)
        df['source_file'] = source_name
        return df

    def generate_source_file_taxon_covered_bar_chart(self, show_legend=False):
        """
        Generate a bar chart for taxon covered bases.
        """
        df_test = self.read_and_append_source(self.test_file_path, 'test file')
        df_control = self.read_and_append_source(self.control_file_path, 'control file')
        df = pd.concat([df_test, df_control])
        df = df[df['taxon_id'] != '*']
        # Use .loc to avoid SettingWithCopyWarning
        df.loc[:, 'taxon_label'] = pd.factorize(df['taxon_id'])[0]

        fig = go.Figure()
        for _, row in df.iterrows():
            fig.add_trace(go.Bar(
                x=[row['source_file']],
                y=[row['taxon_%_covered_bases']],
                name=row['taxon_id'],
                marker=dict(color=self.color_scale[row['taxon_label'] % len(self.color_scale)]),
                hovertemplate='<b>Source File: %{x}</b><br>' +
                              'Estimated Genome Size: %{customdata[1]}<br>' +
                              'Taxon ID: %{text}<br>' +
                              'Taxon Length: %{customdata[0]}<br>' +
                              'Covered Bases: %{y}<br>',
                customdata=[[row['taxon_length'], row['est_genome_size']]],
                text=[row['taxon_id']]
            ))
        fig.update_layout(
            xaxis_title='Source File',
            yaxis_title='Taxon % Covered Bases',
            xaxis=dict(tickmode='array', tickvals=[-0.2, 1.2], ticktext=['test file', 'control file']),
            yaxis=dict(showgrid=True, gridcolor='lightgray'),
            showlegend=show_legend
        )
        self.save_plot_to_html(fig, "source_file_taxon_covered_bar_chart.html")

    def generate_mean_read_length_chart(self):
        """
        Generate a bar chart comparing the taxon mean read length for test and control.
        This plot mimics the taxon covered bases chart, but uses the 'taxon_mean_read_length' column.
        """
        df_test = self.read_and_append_source(self.test_file_path, 'test file')
        df_control = self.read_and_append_source(self.control_file_path, 'control file')
        df = pd.concat([df_test, df_control])
        df = df[df['taxon_id'] != '*']
        df.loc[:, 'taxon_label'] = pd.factorize(df['taxon_id'])[0]

        fig = go.Figure()
        for _, row in df.iterrows():
            fig.add_trace(go.Bar(
                x=[row['source_file']],
                y=[row['taxon_mean_read_length']],
                name=row['taxon_id'],
                marker=dict(color=self.color_scale[row['taxon_label'] % len(self.color_scale)]),
                hovertemplate='<b>Source File: %{x}</b><br>' +
                              'Taxon ID: %{text}<br>' +
                              'Taxon Length: %{customdata[0]}<br>' +
                              'Taxon Mean Read Length: %{y}<br>',
                customdata=[[row['taxon_length']]],
                text=[row['taxon_id']]
            ))
        fig.update_layout(
            xaxis_title='Source File',
            yaxis_title='Taxon Mean Read Length',
            xaxis=dict(tickmode='array', tickvals=[-0.2, 1.2], ticktext=['test file', 'control file']),
            yaxis=dict(showgrid=True, gridcolor='lightgray')
        )
        self.save_plot_to_html(fig, "taxon_mean_read_length_comparison.html")

    def generate_mean_coverage_chart(self):
        """
        Generate a bar chart comparing the taxon mean coverage for test and control.
        This plot mimics the taxon covered bases chart, but uses the 'taxon_mean_coverage' column.
        """
        df_test = self.read_and_append_source(self.test_file_path, 'test file')
        df_control = self.read_and_append_source(self.control_file_path, 'control file')
        df = pd.concat([df_test, df_control])
        df = df[df['taxon_id'] != '*']
        df.loc[:, 'taxon_label'] = pd.factorize(df['taxon_id'])[0]

        fig = go.Figure()
        for _, row in df.iterrows():
            fig.add_trace(go.Bar(
                x=[row['source_file']],
                y=[row['taxon_mean_coverage']],
                name=row['taxon_id'],
                marker=dict(color=self.color_scale[row['taxon_label'] % len(self.color_scale)]),
                hovertemplate='<b>Source File: %{x}</b><br>' +
                              'Taxon ID: %{text}<br>' +
                              'Taxon Length: %{customdata[0]}<br>' +
                              'Taxon Mean Coverage: %{y}<br>',
                customdata=[[row['taxon_length']]],
                text=[row['taxon_id']]
            ))
        fig.update_layout(
            xaxis_title='Source File',
            yaxis_title='Taxon Mean Coverage',
            xaxis=dict(tickmode='array', tickvals=[-0.2, 1.2], ticktext=['test file', 'control file']),
            yaxis=dict(showgrid=True, gridcolor='lightgray')
        )
        self.save_plot_to_html(fig, "taxon_mean_coverage_comparison.html")

    def save_plot_to_html(self, fig, file_name):
        """
        Save a Plotly figure as an HTML file.
        """
        output_file_path = os.path.join(self.output_dir, f"{self.output_prefix}_{file_name}")
        fig.write_html(output_file_path)
        if not self.check_file(output_file_path):
            raise ValueError(f"File {output_file_path} was not created properly.")

    def check_file(self, file_path):
        """
        Check if the file exists and is not empty.
        """
        return os.path.isfile(file_path) and os.path.getsize(file_path) > 0
