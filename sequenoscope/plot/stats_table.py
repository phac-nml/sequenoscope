import pandas as pd
from scipy import stats

class MakeStatsTable:
    def __init__(self, test_file, control_file):
        self.test_data = pd.read_csv(test_file, delimiter='\t')
        self.control_data = pd.read_csv(control_file, delimiter='\t')
        self.parameters = ['est_genome_size', 'est_kmer_coverage_depth', 'total_bases', 'total_fastp_bases',
                            'mean_read_length', 'mean_read_length_reverse', 'taxon_length', 'taxon_%_covered_bases',
                            'taxon_mean_read_length']
        self.result_df = pd.DataFrame(columns=['taxon_id', 'sample_id', 'Parameter', 'Test_Value', 'Control_Value', 
                                               'Ratio', 'Statistical_Test', 'P-Value', 'Significance'])

    def generate_stats(self):
        for index, test_row in self.test_data.iterrows():
            taxon_id = test_row['taxon_id']
            sample_id = test_row['sample_id']

            for param in self.parameters:
                test_value = test_row[param]
                control_value = self.control_data.iloc[index][param]
                ratio = test_value / control_value

                test_normality = stats.shapiro(self.test_data[param])
                control_normality = stats.shapiro(self.control_data[param])

                if test_normality[1] > 0.05 and control_normality[1] > 0.05:
                    t_statistic, p_value = stats.ttest_ind(self.test_data[param], self.control_data[param])
                    statistical_test = 't-test'
                else:
                    u_statistic, p_value = stats.mannwhitneyu(self.test_data[param], self.control_data[param])
                    statistical_test = 'Mann-Whitney U'

                significance = 'Significant' if p_value < 0.05 else 'Not Significant'

                self.result_df = self.result_df.append({
                    'taxon_id': taxon_id,
                    'sample_id': sample_id,
                    'Parameter': param,
                    'Test_Value': test_value,
                    'Control_Value': control_value,
                    'Ratio': ratio,
                    'Statistical_Test': statistical_test,
                    'P-Value': p_value,
                    'Significance': significance,
                }, ignore_index=True)

            self.result_df = self.result_df.append(pd.Series(['-' for _ in self.result_df.columns], index=self.result_df.columns), ignore_index=True)

    def save_to_csv(self, filename='result.csv'):
        self.result_df.to_csv(filename, index=False)
