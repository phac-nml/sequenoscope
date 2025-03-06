#!/usr/bin/env python
import pandas as pd
import plotly.graph_objects as go
import os

class DecisionBarBuilder():
    def __init__(self):
        """Base class for decision bar charts."""
        pass

    def generate_chart(self):
        self.process_data()
        self.create_trace()
        self.create_chart()


class IndependentDecisionStackedBarChart(DecisionBarBuilder):
    def __init__(self, data_path, output_dir, output_prefix="sample", time_bin_unit="seconds"):
        """
        Constructor for the IndependentDecisionStackedBarChart class.
        
        Arguments:
            data_path (str): Path to the data file.
            output_dir (str): Directory where output should be saved.
            output_prefix (str, optional): Prefix for the output file name. Defaults to "sample".
            time_bin_unit (str, optional): Unit of time for binning. Defaults to "seconds".
        """
        self.data_path = data_path
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.classes = {
            "stop_receiving": ["signal_positive"],
            "unblocked": ["data_service_unblock_mux_change"],
            "no_decision": ["signal_negative", "unblock_mux_change"]
        }
        self.time_bin_unit = time_bin_unit
        self.class_order = ["stop_receiving", "unblocked", "no_decision"]

    def process_data(self):
        """
        Load the data, convert start times, build a complete grid of time and decision, 
        and compute the percentage for each decision at each time.
        """
        data = pd.read_csv(self.data_path, sep='\t')
        data['start_time'] = pd.to_datetime(data['start_time'])
        self.total_count_2 = data.groupby('start_time').size().reset_index(name='total_count')
        all_decisions = pd.DataFrame({'decision': ['no_decision', 'stop_receiving', 'unblocked']})
        all_times = pd.DataFrame({'start_time': data['start_time'].unique()})
        complete_grid = all_times.assign(key=1).merge(all_decisions.assign(key=1), on='key').drop('key', axis=1)
        data['decision'] = data['decision'].apply(lambda x: 'no_decision' if x in ['signal_negative', 'unblock_mux_change'] else x)
        decision_count = data.groupby(['start_time', 'decision']).size().reset_index(name='count')
        decision_count['decision'] = decision_count['decision'].map(
            lambda x: next((k for k, v in self.classes.items() if x in v), x))
        decision_count = complete_grid.merge(decision_count, on=['start_time', 'decision'], how='left').fillna({'count': 0})
        total_count = decision_count.groupby('start_time')['count'].sum().reset_index(name='total_count')
        decision_count = decision_count.merge(total_count, on='start_time')
        # Safeguard division: if total_count is 0, percentage becomes 0.
        decision_count['percentage'] = (
            (decision_count['count'] / decision_count['total_count'].replace(0, pd.NA))
            .fillna(0) * 100
        )
        self.decision_count = decision_count

    def create_trace(self):
        """
        Resample the data to compute read counts in the specified time bins.
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
        # Remove any spurious index equal to timestamp 0.
        self.hourly_counts = self.hourly_counts[self.hourly_counts.index != pd.to_datetime(0)]
        self.count_values = self.hourly_counts['read_id'].tolist()

    def convert_time_units(self, x):
        """
        Convert time in seconds to the desired units.
        """
        if self.time_bin_unit in ["5m", "15m", "minutes"]:
            return x / 60
        elif self.time_bin_unit == "hours":
            return x / 3600
        elif self.time_bin_unit == "seconds":
            return x
        else:
            return x

    def create_chart(self):
        """
        Build the stacked bar chart by grouping the decision data into the appropriate time bins
        and overlaying a line trace for read counts.
        """
        self.process_data()
        self.create_trace()
        x_values_all = []

        fig = go.Figure()
        color_palette = ['#2ECC71', '#34495E', '#9B59B6', '#F1C40F']

        for idx, decision in enumerate(self.class_order):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision].copy()
            # Use .loc for assignment
            filtered_data.loc[:, 'start_time_numeric'] = filtered_data['start_time'].astype(int)
            
            # Bin the data according to the time unit
            if self.time_bin_unit == "hours":
                filtered_data.loc[:, 'hourly_bins'] = (filtered_data['start_time_numeric'] // 3600) * 3600
                aggregated_df = filtered_data.groupby('hourly_bins', as_index=False).agg({
                    'count': 'sum', 'total_count': 'sum'})
            elif self.time_bin_unit == "minutes":
                offset = filtered_data['start_time_numeric'].min() % 60
                filtered_data.loc[:, 'minute_bins'] = ((filtered_data['start_time_numeric'] - offset) // 60) * 60
                aggregated_df = filtered_data.groupby('minute_bins', as_index=False).agg({
                    'count': 'sum', 'total_count': 'sum'})
            elif self.time_bin_unit == "5m":
                filtered_data.loc[:, '5m_bins'] = (filtered_data['start_time_numeric'] // 300) * 300
                aggregated_df = filtered_data.groupby('5m_bins', as_index=False).agg({
                    'count': 'sum', 'total_count': 'sum'})
            elif self.time_bin_unit == "15m":
                filtered_data.loc[:, '15m_bins'] = (filtered_data['start_time_numeric'] // 900) * 900
                aggregated_df = filtered_data.groupby('15m_bins', as_index=False).agg({
                    'count': 'sum', 'total_count': 'sum'})
            elif self.time_bin_unit == "seconds":
                aggregated_df = filtered_data[['start_time_numeric', 'percentage']].copy()
                aggregated_df.rename(columns={'start_time_numeric': 'bins', 'percentage': 'percentage'}, inplace=True)
            else:
                aggregated_df = filtered_data[['start_time_numeric', 'percentage']].copy()
                aggregated_df.rename(columns={'start_time_numeric': 'bins', 'percentage': 'percentage'}, inplace=True)

            # Calculate percentage with safeguards
            if self.time_bin_unit in ["hours", "minutes", "5m", "15m"]:
                aggregated_df.loc[:, 'percentage'] = (
                    (aggregated_df['count'] / aggregated_df['total_count'].replace(0, pd.NA))
                    .fillna(0) * 100
                )
                x_data = aggregated_df.iloc[:, 0]  # the binned column (hourly_bins, etc.)
                y_data = aggregated_df['percentage']
            else:
                x_data = aggregated_df['bins']
                y_data = aggregated_df['percentage']

            x_values_all.extend(x_data.map(self.convert_time_units).tolist())
            fig.add_trace(go.Bar(
                x=x_data.map(self.convert_time_units),
                y=y_data,
                name=decision,
                marker_color=color_palette[idx % len(color_palette)],
                yaxis='y'
            ))
        
        x_values_all = list(pd.unique(x_values_all))
        # Add read count trace
        if self.time_bin_unit in ['hours', 'minutes', '5m', '15m']:
            fig.add_trace(go.Scatter(
                x=x_values_all,
                y=self.hourly_counts['read_id'],
                name='Read Count',
                mode='lines',
                line=dict(color='black'),
                visible=True,
                yaxis='y2'
            ))
        elif self.time_bin_unit == 'seconds':
            self.total_count_2.loc[:, 'start_time_numeric'] = self.total_count_2['start_time'].astype(int)
            fig.add_trace(go.Scatter(
                x=self.total_count_2['start_time_numeric'].map(self.convert_time_units),
                y=self.total_count_2['total_count'],
                name='Read Count',
                mode='lines',
                line=dict(color='black'),
                visible=True,
                yaxis='y2'
            ))
        fig.update_layout(
            barmode='stack',
            plot_bgcolor='white',
            xaxis=dict(showgrid=False, title=f'Start Time ({self.time_bin_unit})'),
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
        output_file_path = os.path.join(self.output_dir, f"{self.output_prefix}_independent_decision_bar_chart.html")
        fig.write_html(output_file_path)
        if not os.path.isfile(output_file_path) or os.path.getsize(output_file_path) == 0:
            raise ValueError(f"File {output_file_path} was not created properly.")


class CumulativeDecisionBarChart(DecisionBarBuilder):
    def __init__(self, data_path, output_dir, output_prefix="sample", time_bin_unit="seconds"):
        """
        Constructor for the CumulativeDecisionBarChart class.
        
        Arguments:
            data_path (str): Path to the data file.
            output_dir (str): Directory where output should be saved.
            output_prefix (str, optional): Prefix for the output file name. Defaults to "sample".
            time_bin_unit (str, optional): Unit of time for binning. Defaults to "seconds".
        """
        self.data_path = data_path
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.classes = {
            "stop_receiving": ["signal_positive"],
            "unblocked": ["data_service_unblock_mux_change"],
            "no_decision": ["signal_negative", "unblock_mux_change"]
        }
        self.time_bin_unit = time_bin_unit
        self.class_order = ["stop_receiving", "unblocked", "no_decision"]

    def process_data(self):
        """
        Process the data to compute cumulative counts and percentages per decision over time.
        """
        data = pd.read_csv(self.data_path, sep='\t')
        data['start_time'] = pd.to_datetime(data['start_time'])
        self.total_count_2 = data.groupby('start_time').size().reset_index(name='total_count')
        all_decisions = pd.DataFrame({'decision': ['no_decision', 'stop_receiving', 'unblocked']})
        all_times = pd.DataFrame({'start_time': data['start_time'].unique()})
        complete_grid = all_times.assign(key=1).merge(all_decisions.assign(key=1), on='key').drop('key', axis=1)
        data['decision'] = data['decision'].apply(lambda x: 'no_decision' if x in ['signal_negative', 'unblock_mux_change'] else x)
        decision_count = data.groupby(['start_time', 'decision']).size().reset_index(name='count')
        decision_count['decision'] = decision_count['decision'].map(
            lambda x: next((k for k, v in self.classes.items() if x in v), x))
        decision_count = complete_grid.merge(decision_count, on=['start_time', 'decision'], how='left').fillna({'count': 0})
        decision_count.loc[:, 'cumulative_count'] = decision_count.groupby('decision')['count'].cumsum()
        decision_count.loc[:, 'cumulative_total_count'] = decision_count.groupby('start_time')['cumulative_count'].transform('sum')
        decision_count.loc[:, 'percentage'] = (
            (decision_count['cumulative_count'] / decision_count['cumulative_total_count'].replace(0, pd.NA))
            .fillna(0) * 100
        )
        self.decision_count = decision_count

    def create_trace(self):
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
        if self.time_bin_unit in ["5m", "15m", "minutes"]:
            return x / 60
        elif self.time_bin_unit == "hours":
            return x / 3600
        elif self.time_bin_unit == "seconds":
            return x
        else:
            return x

    def create_chart(self):
        """
        Create a cumulative stacked bar chart with percentage traces and an overlaid line for cumulative read count.
        """
        self.process_data()
        self.create_trace()

        fig = go.Figure()
        color_palette = ['#2ECC71', '#34495E', '#9B59B6', '#F1C40F']

        for idx, decision in enumerate(self.class_order):
            filtered_data = self.decision_count[self.decision_count['decision'] == decision].copy()
            filtered_data.loc[:, 'start_time_numeric'] = filtered_data['start_time'].astype(int)
            
            if self.time_bin_unit == "hours":
                filtered_data.loc[:, 'hourly_bins'] = (filtered_data['start_time_numeric'] // 3600) * 3600
                aggregated_df = filtered_data.groupby('hourly_bins', as_index=False).agg({
                    'cumulative_count': 'sum', 'cumulative_total_count': 'sum'})
            elif self.time_bin_unit == "minutes":
                offset = filtered_data['start_time_numeric'].min() % 60
                filtered_data.loc[:, 'minute_bins'] = ((filtered_data['start_time_numeric'] - offset) // 60) * 60
                aggregated_df = filtered_data.groupby('minute_bins', as_index=False).agg({
                    'cumulative_count': 'sum', 'cumulative_total_count': 'sum'})
            elif self.time_bin_unit == "5m":
                filtered_data.loc[:, '5m_bins'] = (filtered_data['start_time_numeric'] // 300) * 300
                aggregated_df = filtered_data.groupby('5m_bins', as_index=False).agg({
                    'cumulative_count': 'sum', 'cumulative_total_count': 'sum'})
            elif self.time_bin_unit == "15m":
                filtered_data.loc[:, '15m_bins'] = (filtered_data['start_time_numeric'] // 900) * 900
                aggregated_df = filtered_data.groupby('15m_bins', as_index=False).agg({
                    'cumulative_count': 'sum', 'cumulative_total_count': 'sum'})
            elif self.time_bin_unit == "seconds":
                aggregated_df = filtered_data[['start_time_numeric', 'percentage']].copy()
                aggregated_df.rename(columns={'start_time_numeric': 'bins'}, inplace=True)
            else:
                aggregated_df = filtered_data[['start_time_numeric', 'percentage']].copy()
                aggregated_df.rename(columns={'start_time_numeric': 'bins'}, inplace=True)
            
            # Compute percentage with safeguard
            if self.time_bin_unit in ["hours", "minutes", "5m", "15m"]:
                aggregated_df.loc[:, 'percentage'] = (
                    (aggregated_df['cumulative_count'] / aggregated_df['cumulative_total_count'].replace(0, pd.NA))
                    .fillna(0) * 100
                )
                x_data = aggregated_df.iloc[:, 0]
                y_data = aggregated_df['percentage']
            else:
                x_data = aggregated_df['bins']
                y_data = aggregated_df['percentage']

            fig.add_trace(go.Bar(
                x=x_data.map(self.convert_time_units),
                y=y_data,
                name=decision,
                marker_color=color_palette[idx % len(color_palette)],
                yaxis='y'
            ))
        
        self.total_count_2.loc[:, 'start_time_numeric'] = self.total_count_2['start_time'].astype(int)
        cumulative_counts = self.total_count_2['total_count'].cumsum()
        fig.add_trace(go.Scatter(
            x=self.total_count_2['start_time_numeric'].map(self.convert_time_units),
            y=cumulative_counts,
            name='Read Count',
            mode='lines',
            line=dict(color='black'),
            visible=True,
            yaxis='y2'
        ))
        fig.update_layout(
            barmode='stack',
            plot_bgcolor='white',
            xaxis=dict(showgrid=False, title=f'Start Time ({self.time_bin_unit})'),
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
        output_file_path = os.path.join(self.output_dir, f"{self.output_prefix}_cumulative_decision_bar_chart.html")
        fig.write_html(output_file_path)
        if not os.path.isfile(output_file_path) or os.path.getsize(output_file_path) == 0:
            raise ValueError(f"File {output_file_path} was not created properly.")

    def check_files(self, files_to_check):
        if isinstance(files_to_check, str):
            files_to_check = [files_to_check]
        for f in files_to_check:
            if not os.path.isfile(f) or os.path.getsize(f) == 0:
                return False
        return True
