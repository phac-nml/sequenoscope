#!/usr/bin/env python
import os
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go

class SeqManifestPlotter:
    # Class level attributes
    test_file_path = None
    control_file_path = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None
    
    def __init__(self, test_file_path, control_file_path, output_dir, output_prefix="sample"):
        """
        Initialize the SeqManifestPlotter with file paths and output details.
        
        Arguments:
            test_file_path: str
                Path to the test file
            control_file_path: str
                Path to the control file
            output_dir: str
                Directory path where the outputs will be saved
            output_prefix: str, optional
                Prefix for the output files (default is "sample")
        """
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.test_file_path = test_file_path
        self.control_file_path = control_file_path
        self.color_scale = [
        '#FF7F0E',  # Original
        '#1F77B4',  # Original
        '#FFC0CB',  # Original
        '#2CA02C',  # Original
        '#D62728',  # Original
        '#9467BD',  # Original
        '#8C564B',  # Brown
        '#E377C2',  # Pink
        '#7F7F7F',  # Gray
        '#BCBD22',  # Lime
        '#17BECF',  # Cyan
        '#1A55FF',  # Blue
        '#FF991A',  # Orange
        '#B2912F',  # Bronze
        '#197319',  # Dark Green
        '#DB2727',  # Red
        '#9C27B0',  # Purple
        '#03A9F4',  # Light Blue
        '#FFEB3B',  # Yellow
        '#009688'   # Teal
        ]

    def read_data_csv(self, path):
        """
        Reads a CSV file and returns its content as a DataFrame.
        
        Arguments:
            path: str
                Path to the CSV file

        Returns:
            pd.DataFrame:
                Data from the file as a pandas DataFrame
        """
        try:
            df = pd.read_csv(path, delimiter='\t')
        except Exception as e:
            print(f"Error reading data from file: {path}. Error: {e}")
            df = pd.DataFrame()
        return df

    def read_and_append_source(self, path, source_name):
        """
        Reads data from a CSV file and appends a source name to the DataFrame.

        Arguments:
            path: str
                Path to the CSV file
            source_name: str
                Name to append as source_file

        Returns:
            pd.DataFrame:
                DataFrame with appended source_file column
        """
        df = self.read_data_csv(path)
        df['source_file'] = source_name
        return df

    def generate_source_file_taxon_covered_bar_chart(self, show_legend=False):
        """
        Generate a bar chart that visualizes taxon covered bases for test and control files.
        """
        df1 = self.read_and_append_source(self.test_file_path, 'test file')
        df2 = self.read_and_append_source(self.control_file_path, 'control file')
        df = pd.concat([df1, df2])

        df = df[df['taxon_id'] != '*']
        df['taxon_label'] = pd.factorize(df['taxon_id'])[0]
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
            xaxis=dict(tickmode='array', tickvals=[-0.2, 1.2], ticktext=['test_file', 'control_file']),
            yaxis=dict(showgrid=True, gridcolor='lightgray'),
            showlegend=show_legend
        )
        self.save_plot_to_html(fig, "source_file_taxon_covered_bar_chart.html")

    def generate_box_plot(self, parameter):
        """
        Generate a box plot comparing a specific parameter from test and control files.
        
        Arguments:
            parameter: str
                Parameter/column name from the data to be visualized
        """
        df1 = self.read_data_csv(self.test_file_path)
        df2 = self.read_data_csv(self.control_file_path)
        fig = go.Figure()
        fig.add_trace(go.Box(
            y=df1[parameter],
            name='Test Sample',
            showlegend=False 
        ))
        fig.add_trace(go.Box(
            y=df2[parameter],
            name='Control Sample',
            showlegend=False 
        ))
        fig.update_layout(yaxis_title=parameter)
        self.save_plot_to_html(fig, "box_plot.html")

    def generate_ratio_bar_chart(self):
        test_data = self.read_data_csv(self.test_file_path)
        control_data = self.read_data_csv(self.control_file_path)

        if not test_data.columns.equals(control_data.columns):
            print("ERROR: columns for test and control are not equivalent. Please Re-run analyze module with the same parameters.")
            sys.exit()
        
        # Exclude unwanted columns for calculation
        parameter_names = test_data.columns.difference(['sample_id', 'taxon_id'])
        
        # Calculate the average for each parameter
        avg_test_values = test_data[parameter_names].mean()
        avg_control_values = control_data[parameter_names].mean()
        
        # Calculate the ratio of averages, handling zero in denominator
        ratio_values = avg_test_values / avg_control_values.where(avg_control_values != 0, np.nan)
        
        fig = go.Figure()
        for i, parameter in enumerate(parameter_names):
            fig.add_trace(go.Bar(
                x=[parameter],
                y=[ratio_values[parameter]],
                name='Ratio',
                marker=dict(color=self.color_scale[i % len(self.color_scale)]),
                hovertemplate='<b>Parameter:</b> %{x}<br>' +
                            '<b>Ratio:</b> %{y:.2f}<br>' +
                            '<b>Average Test Value:</b> %{customdata[0]:.2f}<br>' +
                            '<b>Average Control Value:</b> %{customdata[1]:.2f}',
                customdata=[[avg_test_values[parameter], avg_control_values[parameter]]],
                showlegend=False
            ))
        fig.update_layout(
            xaxis_title='Parameter',
            yaxis_title='Ratio (Average Test / Average Control)',
            font_family='Arial',
            font_color='rgb(64, 64, 64)',
            margin=dict(l=80, r=80, t=40, b=50),
            height=500
        )
        self.save_plot_to_html(fig, "ratio_bar_chart.html")

    def generate_single_ratio_bar_chart(self, parameter):
        test_data = self.read_data_csv(self.test_file_path)
        control_data = self.read_data_csv(self.control_file_path)
        
        # Calculate the average for the specified parameter
        avg_test_value = test_data[parameter].mean()
        avg_control_value = control_data[parameter].mean()
        
        # Calculate the ratio of averages, handling zero in denominator
        ratio_value = avg_test_value / avg_control_value if avg_control_value != 0 else np.nan
        
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=[parameter],
            y=[ratio_value],
            name='Ratio',
            marker_color='#1f77b4',
            hovertemplate='<b>Parameter:</b> %{x}<br><b>Ratio:</b> %{y:.2f}',
        ))
        fig.update_layout(
            xaxis_title='Parameter',
            yaxis_title='Ratio (Average Test / Average Control)',
            title={
                'text': 'Ratio of Average Test to Average Control for ' + parameter,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            showlegend=False,
            plot_bgcolor='rgba(0,0,0,0)',
            font_family='Arial',
            font_color='rgb(64, 64, 64)',
            title_font=dict(size=24),
            margin=dict(l=50, r=50, t=80, b=50),
            height=300
        )
        self.save_plot_to_html(fig, "single_ratio_bar_chart.html")

    def save_plot_to_html(self, fig, file_name):
        """
        Save a plotly figure as an HTML file in the specified output directory.

        Arguments:
            fig: plotly.graph_objects.Figure
                Figure to save
            file_name: str
                Name for the output file
        """
        output_file_path = os.path.join(self.output_dir, self.output_prefix + "_" + file_name)
        fig.write_html(output_file_path)

        self.status = self.check_files(output_file_path)
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

    def check_files(self, files_to_check):
        """
        check if the output file exists and is not empty

        Arguments:
            files_to_check: list
                list of file paths

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        if isinstance (files_to_check, str):
            files_to_check = [files_to_check]
        for f in files_to_check:
            if not os.path.isfile(f):
                return False
            elif os.path.getsize(f) == 0:
                return False
        return True 

