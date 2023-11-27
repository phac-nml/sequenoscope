#!/usr/bin/env python
import pandas as pd
import os
from scipy import stats

class MakeStatsTable:
    # Class attributes to store various data paths, statuses, and error messages.
    test_file = None
    control_file = None
    status = False
    error_messages = None
    output_dir = None
    output_prefix = None

    def __init__(self, test_file, control_file, output_dir, output_prefix="sample"):
        """
        Constructor for the MakeStatsTable class.

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

        self.parameters = self.test_data.columns.difference(['sample_id', 'taxon_id', 'taxon_length'])

        # self.parameters = ['est_genome_size', 'est_coverage', 'total_bases', 'total_fastp_bases',
        #                    'mean_read_length', 'taxon_length',
        #                    'taxon_%_covered_bases', 'taxon_mean_read_length']
        
        self.result_df = pd.DataFrame(columns=['Parameter', 'Test_Value', 'Control_Value', 'Ratio',
                                               'Statistical_Test', 'P-Value', 'Significance', 'taxon_id'])

    def generate_stats(self):
        for index, test_row in self.test_data.iterrows():
            taxon_id = test_row['taxon_id']
            sample_id = test_row['sample_id']
            for param in self.parameters:
                test_value = test_row[param]
                control_value = self.control_data.iloc[index][param]
                ratio = test_value / control_value

                try:
                    test_normality = stats.shapiro(self.test_data[param])
                    control_normality = stats.shapiro(self.control_data[param])
                except ValueError:
                    # Setting the values as undefined in case of an error
                    test_normality = ('undefined', 'N/A')
                    control_normality = ('undefined', 'N/A')

                # If undefined, skip statistical tests
                if test_normality[1] == 'N/A' or control_normality[1] == 'N/A':
                    statistical_test = 'undefined'
                    p_value = 'N/A'
                    significance = 'undefined'
                else:
                    if test_normality[1] > 0.05 and control_normality[1] > 0.05:
                        t_statistic, p_value = stats.ttest_ind(self.test_data[param], self.control_data[param])
                        statistical_test = 't-test'
                    else:
                        u_statistic, p_value = stats.mannwhitneyu(self.test_data[param], self.control_data[param])
                        statistical_test = 'Mann-Whitney U'
                    significance = 'Significant' if p_value < 0.05 else 'Not Significant'

                new_row = pd.DataFrame([{'taxon_id': taxon_id,
                                        'sample_id': sample_id,
                                        'Parameter': param,
                                        'Test_Value': test_value,
                                        'Control_Value': control_value,
                                        'Ratio': ratio,
                                        'Statistical_Test': statistical_test,
                                        'P-Value': p_value,
                                        'Significance': significance}])

                self.result_df = pd.concat([self.result_df, new_row], ignore_index=True)


    def save_to_csv(self, filename='stat_results.csv'):
        """
        Saves the result data frame to a CSV file.

        Arguments:
        - filename: str
            Name of the output file. Default is "stat_results.csv".
        """
        output_file_path = os.path.join(self.output_dir, self.output_prefix + "_" + filename)
        self.result_df.to_csv(output_file_path, index=False)

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