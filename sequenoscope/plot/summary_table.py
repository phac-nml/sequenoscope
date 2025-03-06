#!/usr/bin/env python
import pandas as pd
import os

class SummaryTable:
    def __init__(self, test_file, control_file, output_dir, output_prefix="sample"):
        """
        Constructor for the SummaryTable class.

        Arguments:
          - test_file: str
              Path to the test data file.
          - control_file: str
              Path to the control data file.
          - output_dir: str
              Directory where the output will be saved.
          - output_prefix: str
              Prefix for the output filename. Default is "sample".
        """
        self.test_file = test_file
        self.control_file = control_file
        self.output_dir = output_dir
        self.output_prefix = output_prefix

        self.test_data = pd.read_csv(test_file, delimiter='\t')
        self.control_data = pd.read_csv(control_file, delimiter='\t')

        # Exclude non-parameter columns
        self.parameters = self.test_data.columns.difference(['sample_id', 'taxon_id', 'taxon_length'])

        # This will hold our summary table
        self.summary_df = None

    def generate_summary(self):
        """
        Generates a summary table with columns:
          Parameter, Test_Value, Control_Value, taxon_id
        Each row corresponds to a parameter for a given taxon.
        """
        rows = []
        # Assume that rows in test_data and control_data correspond by index.
        for i, test_row in self.test_data.iterrows():
            taxon_id = test_row['taxon_id']
            for param in self.parameters:
                test_value = test_row[param]
                control_value = self.control_data.iloc[i][param]
                rows.append({
                    'taxon_id': taxon_id,
                    'Parameter': param,
                    'Test_Value': test_value,
                    'Control_Value': control_value
                })
        self.summary_df = pd.DataFrame(rows)

    def save_to_csv(self, filename='summary_table.csv'):
        """
        Saves the summary table to a CSV file.

        Arguments:
          - filename: str
              Name of the output file. Default is "summary_table.csv".
        """
        if self.summary_df is None:
            self.generate_summary()
        output_file_path = os.path.join(self.output_dir, f"{self.output_prefix}_{filename}")
        self.summary_df.to_csv(output_file_path, index=False)

        if not os.path.isfile(output_file_path) or os.path.getsize(output_file_path) == 0:
            raise ValueError(f"File {output_file_path} was not created properly.")
