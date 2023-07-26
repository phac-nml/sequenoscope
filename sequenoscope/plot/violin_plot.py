import pandas as pd
import plotly.express as px

class ViolinPlotter:
    def __init__(self, test_file, control_file, quality_metric='read_quality', fraction=0.1):
        self.test_file = test_file
        self.control_file = control_file
        self.fraction = fraction
        self.quality_metric = quality_metric
        self.data = None

    def process_file(self, file_path, source_file):
        # Read the file and add source_file column
        df = pd.read_csv(file_path, delimiter='\t')
        df['source_file'] = source_file

        # Check if quality_metric column exists in the dataframe
        if self.quality_metric not in df.columns:
            raise ValueError(f"Quality metric '{self.quality_metric}' not found in the dataframe columns.")

        # Find the min and max values
        max_value = df[self.quality_metric].max()
        min_value = df[self.quality_metric].min()

        # Create a DataFrame with min and max rows
        max_row = pd.DataFrame({self.quality_metric: [max_value], 'source_file': [source_file]})
        min_row = pd.DataFrame({self.quality_metric: [min_value], 'source_file': [source_file]})

        # Sample a fraction of the data
        sampled_df = df.groupby('source_file').apply(lambda x: x.sample(frac=self.fraction)).reset_index(drop=True)

        # Append the min and max rows and the sampled rows
        processed_chunks = pd.concat([min_row, max_row, sampled_df])

        return processed_chunks

    def process_files(self):
        # Process the test file
        processed_test = self.process_file(self.test_file, 'Test')

        # Process the control file
        processed_control = self.process_file(self.control_file, 'Control')

        # Combine the chunks into a single DataFrame
        self.data = pd.concat([processed_test, processed_control])

    def create_violin_plot(self):
        if self.data is None:
            raise ValueError("No data to plot. Please run the process_files() method first.")

        # Create the violin plot
        fig = px.violin(self.data, x='source_file', y=self.quality_metric, points=False)

        # Set the title and labels
        fig.update_layout(
            xaxis_title='Sequencing Type',
            yaxis_title='Read Quality Score'
        )

        # Show the plot
        fig.write_html("violin_comparison_plot.html")
