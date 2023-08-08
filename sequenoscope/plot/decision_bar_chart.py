#!/usr/bin/env python
import pandas as pd
import plotly.graph_objects as go
import os

class IndependentDecisionStackedBarChart():
    data_path = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, data_path, output_dir, output_prefix="sample", time_bin_unit="seconds"):
        self.data_path = data_path
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        #make the classes a constat later on
        self.classes = {
            "stop_receiving": ["signal_positive"],
            "unblocked": ["data_service_unblock_mux_change"],
            "no_decision": ["signal_negative", "unblock_mux_change"]
        }
        self.time_bin_unit = time_bin_unit

    def process_data(self):
        data = pd.read_csv(self.data_path, sep='\t')
        data['start_time'] = pd.to_datetime(data['start_time'])
        self.total_count_2 = data.groupby('start_time').size().reset_index(name='total_count')
        decision_count = data.groupby(['start_time', 'decision']).size().reset_index(name='count')
        decision_count['decision'] = decision_count['decision'].apply(lambda x: 'no_decision' if x in ['signal_negative', 'unblock_mux_change'] else x)
        decision_count = decision_count.groupby(['start_time', 'decision']).sum().reset_index()
        decision_count['decision'] = decision_count['decision'].map(lambda x: next((k for k, v in self.classes.items() if x in v), x))
        total_count = decision_count.groupby('start_time')['count'].sum().reset_index(name='total_count')
        decision_count = decision_count.merge(total_count, on='start_time')
        decision_count['percentage'] = decision_count['count'] / decision_count['total_count'] * 100
        self.decision_count = decision_count

    def create_trace(self):
        df = pd.read_csv(self.data_path, sep='\t')
        df['start_time'] = pd.to_datetime(df['start_time'], unit='s')
        df.set_index('start_time', inplace=True)
        if self.time_bin_unit == "hours":
            self.hourly_counts = df.resample('30T').count()
        elif self.time_bin_unit == "minutes":
            self.hourly_counts = df.resample('1T').count()
        elif self.time_bin_unit == "seconds":
            self.hourly_counts = df.resample('1S').count()
        self.hourly_counts = self.hourly_counts[self.hourly_counts.index != pd.to_datetime(0)]
        self.count_values = self.hourly_counts['read_id'].tolist()

    def create_chart(self):
        self.process_data()
        self.create_trace()
        convert_time_units = lambda x: pd.to_timedelta(x, unit='s') / pd.Timedelta('1{}'.format(self.time_bin_unit))

        x_values = []
        fig = go.Figure()

        color_palette = ['blue', 'green', 'orange', 'red'] 
        for idx, decision in enumerate(self.decision_count['decision'].unique()):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)
            if self.time_bin_unit == "hours":
                filtered_data = filtered_data[
                    (filtered_data['start_time_numeric'] % 3600 == 0) | (filtered_data['start_time_numeric'] % 1800 == 0) | (filtered_data['start_time_numeric'] == 0)
                    ]
            elif self.time_bin_unit == "minutes":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 60 == 0)]
            x_values.extend(filtered_data['start_time_numeric'].map(convert_time_units))  
            fig.add_trace(go.Bar(x=filtered_data['start_time_numeric'].map(convert_time_units),
                                 y=filtered_data['percentage'],
                                 name=decision,
                                 marker_color=color_palette[idx],
                                 yaxis='y'))
        x_values= list(pd.unique(x_values))
        if self.time_bin_unit == 'hours' or self.time_bin_unit == 'minutes':
            fig.add_trace(go.Scatter(x=x_values,
                                    y=self.hourly_counts['read_id'],
                                    name='Read Count',
                                    mode='lines',
                                    line=dict(color='black'),
                                    visible=True,
                                    yaxis='y2'))
        elif self.time_bin_unit == 'seconds':
            self.total_count_2['start_time_numeric'] = self.total_count_2['start_time'].astype(int)
            fig.add_trace(go.Scatter(x=self.total_count_2['start_time_numeric'].astype(int).map(convert_time_units),
                                y=self.total_count_2['total_count'],
                                name='Read Count',
                                mode='lines',
                                line=dict(color='black'),
                                visible=True,
                                yaxis='y2'))
        fig.update_layout(
            barmode='stack',
            plot_bgcolor='white',
            xaxis=dict(showgrid=False, title='Start Time ({})'.format(self.time_bin_unit)),
            yaxis=dict(showgrid=False, title='Percentage'),
            yaxis2=dict(showgrid=False, title='Read Count', overlaying='y', side='right'),  
            legend=dict(
                orientation='h',
                yanchor='top',
                y=1.05,
                xanchor='right',
                x=1
            ),
            margin=dict(t=50, b=50)    
        )
        output_file_path = os.path.join(self.output_dir, self.output_prefix + "_independent_decision_bar_chart.html")
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

class CumulativeDecisionBarChart:
    data_path = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, data_path, output_dir, time_bin_unit, output_prefix="sample"):
        self.data_path = data_path
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        #make the classes a constat later on
        self.classes = {
            "stop_receiving": ["signal_positive"],
            "unblocked": ["data_service_unblock_mux_change"],
            "no_decision": ["signal_negative", "unblock_mux_change"]
        }
        self.time_bin_unit = time_bin_unit

    def process_data(self):
        data = pd.read_csv(self.data_path, sep='\t')
        data['start_time'] = pd.to_datetime(data['start_time'])
        self.total_count_2 = data.groupby('start_time').size().reset_index(name='total_count')
        decision_count = data.groupby(['start_time', 'decision']).size().reset_index(name='count')
        decision_count['decision'] = decision_count['decision'].apply(lambda x: 'no_decision' if x in ['signal_negative', 'unblock_mux_change'] else x)
        decision_count = decision_count.groupby(['start_time', 'decision']).sum().reset_index()
        decision_count['decision'] = decision_count['decision'].map(lambda x: next((k for k, v in self.classes.items() if x in v), x))
        decision_count['cumulative_count'] = decision_count.groupby('decision')['count'].cumsum()
        decision_count['cumulative_total_count'] = decision_count.groupby('start_time')['cumulative_count'].transform('sum')
        decision_count['percentage'] = decision_count['cumulative_count'] / decision_count['cumulative_total_count'] * 100
        self.decision_count = decision_count
        self.x_values = []
    
    def create_trace(self):
        df = pd.read_csv(self.data_path, sep='\t')
        df['start_time'] = pd.to_datetime(df['start_time'])
        df.set_index('start_time', inplace=True)
        if self.time_bin_unit == "hours":
            self.hourly_counts = df.resample('30T').count()
        elif self.time_bin_unit == "minutes":
            self.hourly_counts = df.resample('1T').count()
        elif self.time_bin_unit == "seconds":
            self.hourly_counts = df.resample('1S').count()

        self.hourly_counts = self.hourly_counts[self.hourly_counts.index != pd.to_datetime(0)]
        self.count_values = self.hourly_counts['read_id'].tolist()

    def create_chart(self):
        self.process_data()
        self.create_trace()

        convert_time_units = lambda x: pd.to_timedelta(x, unit='s') / pd.Timedelta('1{}'.format(self.time_bin_unit))

        fig = go.Figure()

        # Define a custom color palette for the decisions
        color_palette = ['blue', 'green', 'orange', 'red']  
        for idx, decision in enumerate(self.decision_count['decision'].unique()):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)

            if self.time_bin_unit == "hours":
                filtered_data = filtered_data[
                    (filtered_data['start_time_numeric'] % 3600 == 0) | (filtered_data['start_time_numeric'] % 1800 == 0) | (filtered_data['start_time_numeric'] == 0)
                    ]
            elif self.time_bin_unit == "minutes":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 60 == 0)]
            
            self.x_values.extend(filtered_data['start_time_numeric'].map(convert_time_units))

            fig.add_trace(go.Bar(x=filtered_data['start_time_numeric'].map(convert_time_units),
                                 y=filtered_data['percentage'],
                                 name=decision,
                                 marker_color=color_palette[idx],
                                 yaxis='y'))
        
        self.x_values = list(pd.unique(self.x_values))

        self.total_count_2['start_time_numeric'] = self.total_count_2['start_time'].astype(int)
        cumulative_counts = self.total_count_2['total_count'].cumsum()
        fig.add_trace(go.Scatter(x=self.total_count_2['start_time_numeric'].astype(int).map(convert_time_units),
                                    y=cumulative_counts,
                                    name='Read Count',
                                    mode='lines',
                                    line=dict(color='black'),
                                    visible=True,
                                    yaxis='y2'))

        fig.update_layout(
            barmode='stack',
            plot_bgcolor='white',
            xaxis=dict(showgrid=False, title='Start Time ({})'.format(self.time_bin_unit)),
            yaxis=dict(showgrid=False, title='Percentage'),
            yaxis2=dict(showgrid=False, title='Cumulative Read Count', overlaying='y', side='right'),
            legend=dict(
                orientation='h',
                yanchor='top',
                y=1.05,
                xanchor='right',
                x=1
            ),
            margin=dict(t=50, b=50)
        )
        output_file_path = os.path.join(self.output_dir, self.output_prefix + "_cumulative_decision_bar_chart.html")
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
