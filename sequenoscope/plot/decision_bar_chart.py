#!/usr/bin/env python
import pandas as pd
import plotly.graph_objects as go
import os

class DecisionBarBuilder():
    def __init__(self):
        """
        Constructor for the DecisionBarBuilder class.
        """
        pass

    def generate_chart(self):
        """
        Generates a chart by processing data, creating a trace, and then creating the chart.
        """
        self.process_data()
        self.create_trace()
        self.create_chart()
    

class IndependentDecisionStackedBarChart(DecisionBarBuilder):
    # Class attributes to store various data paths, statuses, and error messages.
    data_path = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, data_path, output_dir, output_prefix="sample", time_bin_unit="seconds"):
        """
        Constructor for the IndependentDecisionStackedBarChart class.
        
        Arguments:
            data_path (str): 
                Path to the data file.
            output_dir (str): 
                Directory where output should be saved.
            output_prefix (str, optional): 
                Prefix for the output file name. Defaults to "sample".
            time_bin_unit (str, optional): 
                Unit of time for binning. Defaults to "seconds".
        """
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
        self.class_order = ["stop_receiving", "unblocked", "no_decision"]

    def process_data(self):
        """
        Process data by loading it, parsing the date-times, and performing various aggregations.
        """
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
        """
        Generate trace based on the data path and the time bin unit provided in the instance.
        """
        df = pd.read_csv(self.data_path, sep='\t')
        df['start_time'] = pd.to_datetime(df['start_time'], unit='s')
        df.set_index('start_time', inplace=True)
        if self.time_bin_unit == "hours":
            self.hourly_counts = df.resample('30T').count()
        elif self.time_bin_unit == "15m":
            self.hourly_counts = df.resample('15T').count()
        elif self.time_bin_unit == "5m":
            self.hourly_counts = df.resample('5T').count()
        elif self.time_bin_unit == "minutes":
            self.hourly_counts = df.resample('1T').count()
        elif self.time_bin_unit == "seconds":
            self.hourly_counts = df.resample('1S').count()
        self.hourly_counts = self.hourly_counts[self.hourly_counts.index != pd.to_datetime(0)]
        self.count_values = self.hourly_counts['read_id'].tolist()

    def convert_time_units(self, x):
        """
        Convert the time units based on the specified bin unit.

        Arguments:
            x (int or float): 
                Numeric value representing time.

        Returns:
            int or float: 
                Converted time value.
        """
        if self.time_bin_unit == "5m" or self.time_bin_unit == "15m":
            return x / 60  # convert seconds to minutes
        elif self.time_bin_unit == "hours":
            return x / 3600  # convert seconds to hours
        elif self.time_bin_unit == "minutes":
            return x / 60
        elif self.time_bin_unit == "seconds":
            return x
        else:
            return x
        
    def generate_x_values(self):
        """
        Generate x values for the chart based on the decision count.
        """
        self.x_values = []
        
        for decision in self.decision_count['decision'].unique():
            filtered_data = self.decision_count[self.decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)
            if self.time_bin_unit == "hours":
                filtered_data = filtered_data[
                    (filtered_data['start_time_numeric'] % 3600 == 0) | 
                    (filtered_data['start_time_numeric'] % 1800 == 0) | 
                    (filtered_data['start_time_numeric'] == 0)
                ]
            elif self.time_bin_unit == "15m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 900 == 0)]
            elif self.time_bin_unit == "5m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 300 == 0)]
            elif self.time_bin_unit == "minutes":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 60 == 0)]

            self.x_values.extend(filtered_data['start_time_numeric'].map(self.convert_time_units))
        
        self.x_values = list(pd.unique(self.x_values))
        self.x_values.sort()


    def create_chart(self):
        """
        Create a stacked bar chart based on the processed data and x values.
        """
        self.process_data()
        self.create_trace()
        self.generate_x_values()
    
        fig = go.Figure()

        color_palette = ['#2ECC71', '#34495E', '#9B59B6', '#F1C40F']
        for idx, decision in enumerate(self.class_order):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)
            if self.time_bin_unit == "hours":
                filtered_data = filtered_data[
                    (filtered_data['start_time_numeric'] % 3600 == 0) | (filtered_data['start_time_numeric'] % 1800 == 0) | (filtered_data['start_time_numeric'] == 0)
                    ]
            elif self.time_bin_unit == "15m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 900 == 0)]
            elif self.time_bin_unit == "5m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 300 == 0)]
            elif self.time_bin_unit == "minutes":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 60 == 0)]

            fig.add_trace(go.Bar(x=filtered_data['start_time_numeric'].map(self.convert_time_units),
                                 y=filtered_data['percentage'],
                                 name=decision,
                                 marker_color=color_palette[idx],
                                 yaxis='y'))

        if self.time_bin_unit == 'hours' or self.time_bin_unit == 'minutes' or self.time_bin_unit == '5m' or self.time_bin_unit == '15m':
            fig.add_trace(go.Scatter(x=self.x_values,
                                    y=self.hourly_counts['read_id'],
                                    name='Read Count',
                                    mode='lines',
                                    line=dict(color='black'),
                                    visible=True,
                                    yaxis='y2'))
        elif self.time_bin_unit == 'seconds':
            self.total_count_2['start_time_numeric'] = self.total_count_2['start_time'].astype(int)
            fig.add_trace(go.Scatter(x=self.total_count_2['start_time_numeric'].astype(int).map(self.convert_time_units),
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

class CumulativeDecisionBarChart(DecisionBarBuilder):
    # Class attributes to store various data paths, statuses, and error messages.
    data_path = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, data_path, output_dir, output_prefix="sample", time_bin_unit="seconds"):
        """
        Constructor for the CumulativeDecisionBarChart class.
        
        Arguments:
            data_path (str): 
                Path to the data file.
            output_dir (str): 
                Directory where output should be saved.
            output_prefix (str, optional): 
                Prefix for the output file name. Defaults to "sample".
            time_bin_unit (str, optional): 
                Unit of time for binning. Defaults to "seconds".
        """
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
        self.class_order = ["stop_receiving", "unblocked", "no_decision"]

    def process_data(self):
        """
        Processes the data by loading it, parsing date-times, aggregating, and computing percentages.
        """
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
    
    def create_trace(self):
        """
        Generate trace based on the data path and the time bin unit provided in the instance.
        """
        df = pd.read_csv(self.data_path, sep='\t')
        df['start_time'] = pd.to_datetime(df['start_time'], unit='s')
        df.set_index('start_time', inplace=True)
        if self.time_bin_unit == "hours":
            self.hourly_counts = df.resample('30T').count()
        elif self.time_bin_unit == "15m":
            self.hourly_counts = df.resample('15T').count()
        elif self.time_bin_unit == "5m":
            self.hourly_counts = df.resample('5T').count()
        elif self.time_bin_unit == "minutes":
            self.hourly_counts = df.resample('1T').count()
        elif self.time_bin_unit == "seconds":
            self.hourly_counts = df.resample('1S').count()

        self.hourly_counts = self.hourly_counts[self.hourly_counts.index != pd.to_datetime(0)]
        self.count_values = self.hourly_counts['read_id'].tolist()

    def convert_time_units(self, x):
        """
        Convert the time units based on the specified bin unit.

        Arguments:
            x (int or float): 
                Numeric value representing time.

        Returns:
            int or float: 
                Converted time value.
        """
        if self.time_bin_unit == "5m" or self.time_bin_unit == "15m":
            return x / 60  # convert seconds to minutes
        elif self.time_bin_unit == "hours":
            return x / 3600  # convert seconds to hours
        elif self.time_bin_unit == "minutes":
            return x / 60
        elif self.time_bin_unit == "seconds":
            return x
        else:
            return x

    def create_chart(self):
        """
        Creates a stacked bar chart based on the processed data, x values, and percentages.
        """
        self.process_data()
        self.create_trace()

        fig = go.Figure()

        # Define a custom color palette for the decisions
        color_palette = ['#2ECC71', '#34495E', '#9B59B6', '#F1C40F']  
        for idx, decision in enumerate(self.class_order):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)

            if self.time_bin_unit == "hours":
                filtered_data = filtered_data[
                    (filtered_data['start_time_numeric'] % 3600 == 0) | (filtered_data['start_time_numeric'] % 1800 == 0) | (filtered_data['start_time_numeric'] == 0)
                    ]
            elif self.time_bin_unit == "15m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 900 == 0)]
            elif self.time_bin_unit == "5m":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 300 == 0)]
            elif self.time_bin_unit == "minutes":
                filtered_data = filtered_data[(filtered_data['start_time_numeric'] % 60 == 0)]
            
            fig.add_trace(go.Bar(x=filtered_data['start_time_numeric'].map(self.convert_time_units),
                                 y=filtered_data['percentage'],
                                 name=decision,
                                 marker_color=color_palette[idx],
                                 yaxis='y'))
        

        self.total_count_2['start_time_numeric'] = self.total_count_2['start_time'].astype(int)
        cumulative_counts = self.total_count_2['total_count'].cumsum()
        fig.add_trace(go.Scatter(x=self.total_count_2['start_time_numeric'].astype(int).map(self.convert_time_units),
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
