#!/usr/bin/env python
import subprocess
import os

class MashSketcher:
    def __init__(self, out_directory, out_prefix):
        """
        Initialize the class with out_directory, out_prefix.

        Arguments:
            out_directory: str
                A string to the path where the output files will be stored.
            out_prefix: str
                A designation of what the output files will be named.
        """
        self.results = {}
        self.out_directory = out_directory
        self.out_prefix = out_prefix

    def run_mash_sketch(self, input_files):
        """
        Run the MASH sketch process on the input files. Can handle either single or paired input.

        Arguments:
            input_files: list
                List containing paths to the input file(s).

        Returns:
            dict: Dictionary containing results of the MASH sketch process.
        """
        if len(input_files) == 1:
            self.results = self._sketch_single(input_files[0], self.out_prefix)
        elif len(input_files) == 2:
            self.results = self._sketch_paired(input_files[0], input_files[1])
        else:
            raise ValueError("Input files list must contain either one or two files.")
        return self.results

    def _sketch_single(self, input_file, file_prefix):
        """
        Run MASH sketch on a single input file and extracts the genome size and coverage.

        Arguments:
            input_file: str
                Path to the input file.
            file_prefix: str
                Prefix for the output file.

        Returns:
            dict: Dictionary containing genome size and coverage.
        """
        out_file_hash = os.path.join(self.out_directory, f"{file_prefix}_mash_hash")
        cmd = ["mash", "sketch", "-r", input_file, "-o", out_file_hash, "-k", "27", "-m", "3"]
        if file_prefix != "sample":
            cmd.insert(3, "-C")
            cmd.insert(4, file_prefix)
    
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        mash_results = {}

        for line in result.stderr.split('\n'):
            if "Estimated genome size:" in line:
                genome_size = float(line.split(":")[1].strip())
                genome_size = format(genome_size, '.2f')
                mash_results["Genome Size"] = genome_size
            elif "Estimated coverage:" in line:
                coverage = float(line.split(":")[1].strip())
                coverage = format(coverage, '.2f')
                mash_results["Coverage"] = coverage

        return mash_results

    def _sketch_paired(self, input_file1, input_file2):
        """
        Run MASH sketch on paired input files and computes the average genome size and coverage.
        
        Arguments:
            input_file1: str
                Path to the first input file.
            input_file2: str
                Path to the second input file.

        Returns:
            dict: Dictionary with the average genome size and coverage.
        """
        results1 = self._sketch_single(input_file1, f"{self.out_prefix}_1")
        results2 = self._sketch_single(input_file2, f"{self.out_prefix}_2")
        
        avg_genome_size = (float(results1["Genome Size"]) + float(results2["Genome Size"])) / 2
        avg_coverage = (float(results1["Coverage"]) + float(results2["Coverage"])) / 2
        
        avg_genome_size = format(avg_genome_size, '.2f')
        avg_coverage = format(avg_coverage, '.2f')
        
        return {"Genome Size": avg_genome_size, "Coverage": avg_coverage}