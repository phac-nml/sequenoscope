#!/usr/bin/env python
import pandas as pd
import os
import plotly.express as px

class ViolinBuilder():
    def __init__(self):
        """
        Constructor for the ViolinBuilder class.
        """
        pass

    def generate_chart(self):
        """
        Generates a chart by processing files and then creating the plot.
        """
        self.process_files()
        self.create_violin_plot()


class ViolinPlotter(ViolinBuilder):
    # Class attributes to store various data paths, statuses, and error messages.
    test_file = None
    control_file = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, test_file, control_file, output_dir, output_prefix="sample", quality_metric='read_qscore', fraction=0.1):
        """
        Constructor for the ViolinPlotter class.

        Arguments:
        - test_file: str
            Path to the test data file.
        - control_file: str
            Path to the control data file.
        - output_dir: str
            Directory where the output will be saved.
        - output_prefix: str
            Prefix for the output filename. Default is "sample".
        - quality_metric: str
            Quality metric to be used for the violin plot. Default is 'read_qscore'.
        - fraction: float
            Fraction of data to be sampled for plotting. Default is 0.1.
        """
        self.test_file = test_file
        self.control_file = control_file
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.fraction = fraction
        self.quality_metric = quality_metric
        self.data = None

    def process_file(self, file_path, source_file):
        """
        Process a single file to prepare it for the violin plot.

        Arguments:
        - file_path: str
            Path to the file to process.
        - source_file: str
            Source of the file. It's used as a label in the final plot.

        Returns:
        - pd.DataFrame
            Processed dataframe with the required columns for plotting.
        """
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
        """
        Process both test and control files.
        """
        # Process the test file
        processed_test = self.process_file(self.test_file, 'Test')

        # Process the control file
        processed_control = self.process_file(self.control_file, 'Control')

        # Combine the chunks into a single DataFrame
        self.data = pd.concat([processed_test, processed_control])

    def create_violin_plot(self):
        """
        Generate the violin plot based on processed data.
        """
        if self.data is None:
            raise ValueError("No data to plot. Please run the process_files() method first.")

        # Create the violin plot
        fig = px.violin(self.data, x='source_file', y=self.quality_metric, points=False)

        # Set the title and labels
        fig.update_layout(
            xaxis_title='Source',
            yaxis_title= self.quality_metric
        )

        # Show the plot
        output_file_path = os.path.join(self.output_dir, self.output_prefix + "_violin_comparison_plot.html")
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