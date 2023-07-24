import pandas as pd
import plotly.graph_objects as go

class DataUtil:
    @staticmethod
    def read_data(path):
        try:
            df = pd.read_csv(path, delimiter='\t')
        except Exception as e:
            print(f"Error reading data from file: {path}. Error: {e}")
            df = pd.DataFrame()
        return df

class SourceFileTaxonCoveredBarChart:
    @staticmethod
    def read_and_append_source(path, source_name):
        df = DataUtil.read_data(path)
        df['source_file'] = source_name
        return df

    def __init__(self, file1, file2):
        self.df1 = self.read_and_append_source(file1, 'test file')
        self.df2 = self.read_and_append_source(file2, 'control file')
        self.df = pd.concat([self.df1, self.df2])
        self.color_scale = ['#FF7F0E', '#1F77B4', '#FFC0CB', '#2CA02C', '#D62728', '#9467BD']

    def generate_plot(self):
        self.df = self.df[self.df['taxon_id'] != '*']
        self.df['taxon_label'] = pd.factorize(self.df['taxon_id'])[0]
        fig = go.Figure()
        for _, row in self.df.iterrows():
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
        fig.update_traces(marker=dict(colorscale=self.color_scale))
        self.save_plot(fig, "source_file_taxon_covered_bar_chart.html")

    @staticmethod
    def save_plot(fig, file_name):
        fig.write_html(file_name)

class BoxPlot:
    def __init__(self, file1, file2, parameter):
        self.df1 = DataUtil.read_data(file1)
        self.df2 = DataUtil.read_data(file2)
        self.parameter = parameter

    def generate_plot(self):
        fig = go.Figure()
        fig.add_trace(go.Box(
            y=self.df1[self.parameter],
            name='Test Sample',
            showlegend=False 
        ))
        fig.add_trace(go.Box(
            y=self.df2[self.parameter],
            name='Control Sample',
            showlegend=False 
        ))
        fig.update_layout(yaxis_title=self.parameter)
        SourceFileTaxonCoveredBarChart.save_plot(fig, "box_plot.html")

class RatioBarChart:
    def __init__(self, file1, file2):
        self.test_data = DataUtil.read_data(file1).sort_values(by='sample_id')
        self.control_data = DataUtil.read_data(file2).sort_values(by='sample_id')
        self.color_palette = ['#FF7F0E', '#1F77B4', '#FFC0CB', '#2CA02C', '#D62728', '#9467BD']

    def generate_plot(self):
        parameter_names = self.test_data.columns[1:].difference(['sample_id', 'taxon_id'])
        ratio_values = self.test_data[parameter_names] / self.control_data[parameter_names]
        fig = go.Figure()
        for i, parameter in enumerate(parameter_names):
            fig.add_trace(go.Bar(
                x=[parameter],
                y=[ratio_values.loc[0, parameter]],
                name='Ratio',
                marker=dict(color=self.color_palette[i % len(self.color_palette)]),
                hovertemplate='<b>Parameter:</b> %{x}<br>' +
                              '<b>Ratio:</b> %{y:.2f}<br>' +
                              '<b>Test Value:</b> %{customdata[0]:.2f}<br>' +
                              '<b>Control Value:</b> %{customdata[1]:.2f}',
                customdata=[[self.test_data.loc[0, parameter], self.control_data.loc[0, parameter]]],
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
        SourceFileTaxonCoveredBarChart.save_plot(fig, "ratio_bar_chart.html")

class SingleRatioBarChart:
    def __init__(self, file1, file2, parameter):
        self.test_data = DataUtil.read_data(file1).sort_values(by='sample_id')
        self.control_data = DataUtil.read_data(file2).sort_values(by='sample_id')
        self.parameter = parameter

    def generate_plot(self):
        ratio_values = self.test_data[self.parameter] / self.control_data[self.parameter]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=[self.parameter],
            y=[ratio_values[0]],
            name='Ratio',
            marker_color='#1f77b4',
            hovertemplate='<b>Parameter:</b> %{x}<br><b>Ratio:</b> %{y:.2f}',
        ))
        fig.update_layout(
            xaxis_title='Parameter',
            yaxis_title='Ratio (Test/Control)',
            title={
                'text': 'Ratio of Test to Control for ' + self.parameter,
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
        SourceFileTaxonCoveredBarChart.save_plot(fig, "single_ratio_bar_chart.html")
