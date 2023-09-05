#!/usr/bin/env python
import subprocess
import os

class MashSketcher:
    def __init__(self, out_directory, out_prefix):
        self.results = {}
        self.out_directory = out_directory
        self.out_prefix = out_prefix

    def run_mash_sketch(self, input_files):
        if len(input_files) == 1:
            self.results = self._sketch_single(input_files[0], self.out_prefix)
        elif len(input_files) == 2:
            self.results = self._sketch_paired(input_files[0], input_files[1])
        else:
            raise ValueError("Input files list must contain either one or two files.")
        return self.results

    def _sketch_single(self, input_file, file_prefix):
        out_file_hash = os.path.join(self.out_directory, f"{file_prefix}_mash_hash")
        cmd = ["mash", "sketch", "-r", input_file, "-o", out_file_hash, "-k", "27"]
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
        results1 = self._sketch_single(input_file1, f"{self.out_prefix}_1")
        results2 = self._sketch_single(input_file2, f"{self.out_prefix}_2")
        
        avg_genome_size = (float(results1["Genome Size"]) + float(results2["Genome Size"])) / 2
        avg_coverage = (float(results1["Coverage"]) + float(results2["Coverage"])) / 2
        
        avg_genome_size = format(avg_genome_size, '.2f')
        avg_coverage = format(avg_coverage, '.2f')
        
        return {"Genome Size": avg_genome_size, "Coverage": avg_coverage}