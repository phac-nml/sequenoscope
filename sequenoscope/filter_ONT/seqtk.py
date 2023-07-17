#!/usr/bin/env python
from sequenoscope.utils.__init__ import run_command
import os

class SeqtkRunner:
    csv_file = None
    out_dir = None
    out_prefix = None
    status = False
    error_messages = None
    result_files = {"output_fastq":""}
    
    
    def __init__(self, read_set, csv_file, out_dir, out_prefix):
        """
        Initalize the class with read_set, csv, out_dir, and out_prefix

        Arguments:
            read_set: sequence object
                an object that contains the list of sequence files for analysis
            csv_file: str
                a string to the path where the classified csv files are stored
            out_dir: str
                a string to the path where the output files will be stored
            out_prefix: str
                a designation of what the output files will be named
        """
        self.read_set = read_set
        self.csv_file = csv_file
        self.out_dir = out_dir
        self.out_prefix = out_prefix
    
    def subset_fastq(self):
        """
        generate a fastq file of unblocked reads based on the designation of the reads in the adaptive sampling
        report

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """

        output_fastq = os.path.join(self.out_dir,"{}_subset.fastq".format(self.out_prefix))

        self.result_files["output_fastq"] = output_fastq

        cmd = "seqtk subseq {} {} > {}".format(self.read_set.out_files, self.csv_file, output_fastq)

        (self.stdout, self.stderr) = run_command(cmd)
        self.status = self.check_files([output_fastq])
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