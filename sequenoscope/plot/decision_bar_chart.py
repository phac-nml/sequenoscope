import pandas as pd
import plotly.graph_objects as go

class IndependentDecisionStackedBarChart:
    def __init__(self, data_path, time_bin_unit):
        self.data_path = data_path
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
            title='Decision Distribution vs Read Count',    
        )
        fig.write_html("independent_decision_bar_chart.html")

class CumulativeDecisionBarChart:
    def __init__(self, data_path, time_bin_unit):
        self.data_path = data_path
        self.classes = {
            "stop_receiving": ["signal_positive"],
            "unblocked": ["data_service_unblock_mux_change"],
            "no_decision": ["signal_negative", "unblock_mux_change"]
        }
        self.time_bin_unit = time_bin_unit

    def create_chart(self):
        # Load data from the file
        df = pd.read_csv(self.data_path, sep='\t')

        # Convert the 'start_time' column to datetime format
        df['start_time'] = pd.to_datetime(df['start_time'])

        # Calculate the count for each decision at each start time
        decision_count = df.groupby(['start_time', 'decision']).size().reset_index(name='count')

        # Modify the counting logic to count "signal_negative" and "unblock_mux_change" only once for "no_decision"
        decision_count['decision'] = decision_count['decision'].apply(lambda x: 'no_decision' if x in ['signal_negative', 'unblock_mux_change'] else x)
        decision_count = decision_count.groupby(['start_time', 'decision']).sum().reset_index()

        # Map the decisions based on the dictionary
        decision_count['decision'] = decision_count['decision'].map(lambda x: next((k for k, v in self.classes.items() if x in v), x))

        # Calculate the cumulative count for each decision at each start time
        decision_count['cumulative_count'] = decision_count.groupby('decision')['count'].cumsum()

        # Calculate the cumulative total count for each start time
        decision_count['cumulative_total_count'] = decision_count.groupby('start_time')['cumulative_count'].transform('sum')

        # Calculate the percentage for each decision at each start time
        decision_count['percentage'] = decision_count['cumulative_count'] / decision_count['cumulative_total_count'] * 100

        # Define the time bin unit and selected unit
        convert_time_units = lambda x: pd.to_timedelta(x, unit='s') / pd.Timedelta('1{}'.format(self.time_bin_unit))

        # Create lists to store the x-axis values for the bar chart
        x_values = []

        df['start_time'] = pd.to_datetime(df['start_time'])
        df.set_index('start_time', inplace=True)
        if self.time_bin_unit == "hours":
            hourly_counts = df.resample('30T').count()
        elif self.time_bin_unit == "minutes":
            hourly_counts = df.resample('1T').count()
        elif self.time_bin_unit == "seconds":
            hourly_counts = df.resample('1S').count()

        # Create the stacked bar chart
        fig = go.Figure()

        # Define a custom color palette for the decisions
        color_palette = ['blue', 'green', 'orange', 'red']  
        for idx, decision in enumerate(decision_count['decision'].unique()):
            filtered_data = decision_count[decision_count['decision'] == decision]
            filtered_data['start_time_numeric'] = filtered_data['start_time'].astype(int)

            fig.add_trace(go.Bar(x=filtered_data['start_time_numeric'].map(convert_time_units),
                                 y=filtered_data['percentage'],
                                 name=decision,
                                 marker_color=color_palette[idx],
                                 yaxis='y'))

        if self.time_bin_unit == 'hours' or self.time_bin_unit == 'minutes':
            cumulative_counts = hourly_counts['read_id'].cumsum()
            fig.add_trace(go.Scatter(x=x_values,
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

        # Display the stacked bar chart
        fig.write_html("cumulative_decision_bar_chart.html")

