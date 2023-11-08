#!/usr/bin/env python
import os
import pandas as pd
import sys

class BarcodeStatistics:
    input_csv_file = None
    out_dir = None
    out_prefix = None
    status = False
    error_messages = None
    result_files = {"output_csv_file":""}

    def __init__(self, input_csv_file, out_dir, out_prefix):
        """
        Initialize the BarcodeStatistics class.

        Arguments:
            input_csv_file: str
                Path to the input CSV file containing barcode data.
            out_dir: str
                Directory to save the generated statistics CSV file. Defaults to current directory if not provided.
            out_prefix: str
                Prefix to be added to the output file name.

        """
        self.input_csv_file = input_csv_file
        self.out_dir = out_dir
        self.out_prefix = out_prefix

    def generate_statistics(self):
        """
        Generate statistics based on the barcode data from the input CSV file.

        Returns:
            None, but raises a ValueError if there's an issue generating or saving the statistics.
        """
        
        df = pd.read_csv(self.input_csv_file, delimiter='\t')

        if 'barcode_arrangement' not in df.columns:
            try:
                raise ValueError("The 'barcode_arrangement' column was not found in the sequencing summary. Please check your input file.")
            except:
                print("-"*40)
                print("Error: The 'barcode_arrangement' column was not found in the sequencing summary. Please check your input file.")
                print("-"*40)
                sys.exit()

        columns_of_interest = ["channel", "start_time", "duration", 
                            "sequence_length_template", "mean_qscore_template"]
        
        # Filter out columns that are not in the DataFrame
        available_columns = [col for col in columns_of_interest if col in df.columns]

        all_stats = []
        unique_barcodes = [barcode for barcode in df['barcode_arrangement'].unique() if barcode != 'unclassified']
        sorted_barcodes = sorted(unique_barcodes, key=lambda x: int(x[7:]))
        if 'unclassified' in df['barcode_arrangement'].unique():
            sorted_barcodes.append('unclassified')

        for barcode in sorted_barcodes:
            filtered_df = df[df['barcode_arrangement'] == barcode]
            
            stats_dict = {
                "Barcode": barcode,
                "read_number": len(filtered_df)
            }
            
            for col in available_columns:
                col_name = col.replace("mean_qscore_template", "qscore").replace("sequence_length_template", "sequence_length")
                
                stats_dict[f"{col_name}_min"] = filtered_df[col].min() if col in filtered_df else None
                stats_dict[f"{col_name}_max"] = filtered_df[col].max() if col in filtered_df else None
                
                if col != "channel" and col in filtered_df:
                    stats_dict[f"{col_name}_mean"] = filtered_df[col].mean()
            
            all_stats.append(stats_dict)

        all_stats_df = pd.DataFrame(all_stats)

        output_csv_file = os.path.join(self.out_dir, self.out_prefix + "_barcode_statistics.csv")
        all_stats_df.to_csv(output_csv_file, index=False)
        self.result_files["output_csv_file"] = output_csv_file

        # Check the created file
        self.status = self.check_files([output_csv_file])
        if self.status == False:
            self.error_messages = "One or more files was not created or was empty."
            raise ValueError(str(self.error_messages))
        
    def check_files(self, files_to_check):
        """
        Check if the output file exists and is not empty.

        Arguments:
            files_to_check: list
                List of file paths

        Returns:
            bool:
                Returns True if the generated output file is found and not empty, False otherwise.
        """
        if isinstance(files_to_check, str):
            files_to_check = [files_to_check]
        for f in files_to_check:
            if not os.path.isfile(f):
                return False
            elif os.path.getsize(f) == 0:
                return False
        return True