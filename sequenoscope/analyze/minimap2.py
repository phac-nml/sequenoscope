#!/usr/bin/env python

import os
from sequenoscope.constant import DefaultValues
from sequenoscope.utils.__init__ import run_command


class Minimap2Runner:
    read_set = None
    out_dir = None
    out_prefix = None
    ref_database = None
    threads = 1
    kmer_size = 15
    status = False
    error_messages = None
    result_files =  {"sam_output_file":""}
    paired = False

    def __init__(self, read_set, out_dir, ref_database, out_prefix, threads=1, kmer_size=DefaultValues.minimap2_kmer_size):
        """
        Initalize the class with read_set, out_dir, ref_database, and out_prefix

        Arguments:
            read_set: sequence object
                an object that contains the list of sequence files for analysis
            out_dir: str
                a string to the path where the output files will be stored
            ref_database: str
                a string to the path of reference sequence file
            out_prefix: str
                a designation of what the output files will be named
            threads: int
                an integer representing the number of threads utilized for the operation, default is 1
            kmersize: int
                an integer representing the kmer size utilized for the kat filter method, default is 15
        """
        self.read_set = read_set
        self.out_dir = out_dir
        self.out_prefix = out_prefix
        self.ref_database = ref_database
        self.threads = threads
        self.kmer_size = kmer_size
        self.paired = self.read_set.is_paired

    def run_minimap2(self):
        """
        Run the minimap2 program with the designated paramters selected during the initalization of the class.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        sam_file = os.path.join(self.out_dir,f"{self.out_prefix}.sam")
        
        self.result_files["sam_output_file"] = sam_file

        cmd = ["minimap2", "-ax", "-t", f"{self.threads}", "-k", f"{self.kmer_size}", self.ref_database, 
        self.read_set.out_files, ">", sam_file]
        
        if self.paired:
            cmd.insert(2, "sr")
        else:
            cmd.insert(2, "map-ont")

        cmd_string = " ".join(cmd)

        (self.stdout, self.stderr) = run_command(cmd_string)
        self.status = self.check_files([sam_file])
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
        
    
