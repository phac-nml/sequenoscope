#!/usr/bin/env python
import os
import pandas as pd
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
        self.color_scale = ['#FF7F0E', '#1F77B4', '#FFC0CB', '#2CA02C', '#D62728', '#9467BD']

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

    def generate_source_file_taxon_covered_bar_chart(self):
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
            xaxis=dict(showgrid=True, gridcolor='lightgray'),
            yaxis=dict(showgrid=True, gridcolor='lightgray'),
            showlegend=False
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
        """
        Generate a bar chart that visualizes the ratio between test and control samples.
        """
        test_data = self.read_data_csv(self.test_file_path).sort_values(by='sample_id')
        control_data = self.read_data_csv(self.control_file_path).sort_values(by='sample_id')

        parameter_names = test_data.columns[1:].difference(['sample_id', 'taxon_id'])
        ratio_values = test_data[parameter_names] / control_data[parameter_names]
        fig = go.Figure()
        for i, parameter in enumerate(parameter_names):
            fig.add_trace(go.Bar(
                x=[parameter],
                y=[ratio_values.loc[0, parameter]],
                name='Ratio',
                marker=dict(color=self.color_scale[i % len(self.color_scale)]),
                hovertemplate='<b>Parameter:</b> %{x}<br>' +
                              '<b>Ratio:</b> %{y:.2f}<br>' +
                              '<b>Test Value:</b> %{customdata[0]:.2f}<br>' +
                              '<b>Control Value:</b> %{customdata[1]:.2f}',
                customdata=[[test_data.loc[0, parameter], control_data.loc[0, parameter]]],
                showlegend=False
            ))
        fig.update_layout(
            xaxis_title='Parameter',
            yaxis_title='Ratio (Test/Control)',
            font_family='Arial',
            font_color='rgb(64, 64, 64)',
            margin=dict(l=80, r=80, t=40, b=50),
            height=500
        )
        self.save_plot_to_html(fig, "ratio_bar_chart.html")

    def generate_single_ratio_bar_chart(self, parameter):
        """
        Generate a bar chart visualizing the ratio of a single parameter between test and control samples.
        
        Arguments:
            parameter: str
                Parameter/column name from the data to be visualized
        """
        test_data = self.read_data_csv(self.test_file_path).sort_values(by='sample_id')
        control_data = self.read_data_csv(self.control_file_path).sort_values(by='sample_id')
        ratio_values = test_data[parameter] / control_data[parameter]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=[parameter],
            y=[ratio_values[0]],
            name='Ratio',
            marker_color='#1f77b4',
            hovertemplate='<b>Parameter:</b> %{x}<br><b>Ratio:</b> %{y:.2f}',
        ))
        fig.update_layout(
            xaxis_title='Parameter',
            yaxis_title='Ratio (Test/Control)',
            title={
                'text': 'Ratio of Test to Control for ' + parameter,
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

