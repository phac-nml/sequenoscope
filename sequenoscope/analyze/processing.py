#!/usr/bin/env python
import os
from sequenoscope.utils.__init__ import run_command


class SamBamProcessor:
    file = None
    out_dir = None
    out_prefix = None
    ref_database = None
    threads = 1
    status = False
    error_messages = None
    result_files = {"bam_output":"", "fastq_output":"", "bam_output":"", "coverage_tsv":""}
    
    
    def __init__(self, file, out_dir, ref_database, out_prefix, thread=1):
        """
        Initalize the class with read_set, out_dir, ref_database, and out_prefix

        Arguments:
            file: str
                path to sam, bam, or fastq file for processing
            out_dir: str
                a string to the path where the output files will be stored
            ref_database: str
                a string to the path of reference sequence file
            out_prefix: str
                a designation of what the output files will be named
            threads: int
                an integer representing the number of threads utilized for the operation, default is 1
        """
        self.file = file
        self.out_dir = out_dir
        self.ref_database = ref_database
        self.out_prefix = out_prefix
        self.threads = thread

    def run_samtools_bam(self):
        """
        Run the samtools view command on the designated sam file and pipe the output to the samtools sort command for further processing.

        Arguments:
            exclude: bool
                returns all mapped reads if True, and returns all unmapped reads if False

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        bam_output = os.path.join(self.out_dir,f"{self.out_prefix}.bam")
        
        self.result_files["bam_output"] = bam_output
        
        cmd = ["samtools", "view", "-S", "-b", self.file, "|",
               "samtools", "sort", "-@", f"{self.threads}", "-T", self.out_prefix, "--reference", self.ref_database,
               "-o", bam_output]

        cmd_string = " ".join(cmd)

        (self.stdout, self.stderr) = run_command(cmd_string)
        self.status = self.check_files([bam_output])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

    def run_samtools_fastq(self):
        """
        Run the samtools fastq command to convert the designated bam file to a fastq file

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        fastq_output = os.path.join(self.out_dir,f"{self.out_prefix}.fastq")

        self.result_files["fastq_output"] = fastq_output

        cmd = f"samtools fastq {self.file} > {fastq_output}"

        (self.stdout, self.stderr) = run_command(cmd)
        self.status = self.check_files([fastq_output])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))
        
    def run_samtools_import(self):
        """
        Run the samtools fastq command to convert the designated bam file to a fastq file

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        bam_output = os.path.join(self.out_dir,f"{self.out_prefix}.bam")

        self.result_files["bam_output"] = bam_output

        cmd = f"samtools import {self.file} > {bam_output}"

        (self.stdout, self.stderr) = run_command(cmd)
        self.status = self.check_files([bam_output])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

    def run_bedtools(self, nonzero=False):
        """
        Run the bedtools genomcov command on the designated bam file and return a tsv of the depth per base

        Arguments:
            nonzero: bool
                returns only depths greater than zero if True, and returns all depths if False

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        coverage_tsv = os.path.join(self.out_dir, f"{self.out_prefix}.tsv")

        self.result_files["coverage_tsv"] = coverage_tsv

        cmd = ["bedtools", "genomecov", "-ibam", self.file, ">", coverage_tsv]

        if nonzero == True:
            cmd.insert(4, "-dz")
        else:
            cmd.insert(4, "-d")

        cmd_string = " ".join(cmd)

        (self.stdout, self.stderr) = run_command(cmd_string)
        self.status = self.check_files([coverage_tsv])
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
    