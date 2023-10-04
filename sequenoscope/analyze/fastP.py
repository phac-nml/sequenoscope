#!/usr/bin/env python
from sequenoscope.utils.__init__ import run_command
from sequenoscope.utils.sequence_class import Sequence
import os

class FastPRunner:
    read_set = None
    out_dir = None
    out_prefix = None
    min_read_len = 0
    max_read_len = 0
    trim_front_bp = 0
    trim_tail_bp = 0
    report_only = True
    dedup = False
    threads = 1
    status = False
    error_messages = None
    result_files = {"html":"", "json":"", "output_files_fastp":[]}
    paired = False

    def __init__(self, read_set, out_dir, out_prefix, qualified_quality_phred=15, min_read_len=15, max_read_len=0, 
    trim_front_bp=0, trim_tail_bp=0, report_only=True, dedup=True, threads=1):
        """
        Initalize the class with read_set, out_dir, and out_prefix

        Arguments:
            read_set: sequence object
                an object that contains the list of sequence files for analysis
            out_dir: str
                a string to the path where the output files will be stored
            out_prefix: str
                a designation of what the output files will be named
            out_prefix_2: str
                a designation of what the output files will be named in the event that paired end reads are provided
            qualified_quality_phred
                the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
            min_read_len: int
                reads shorter than the integer specified required will be discarded, default is 15
            max_read_len: int
                reads longer than the integer specified required will be discarded, default is 0 meaning no limitation
            trim_front_bp: int
                trimming how many bases in front for read_set, default is 0
            trim_tail_bp: int
                trimming how many bases in tail for read_set, default is 0
            report_only: bool
                a designation of wheather or not to output the fastq file after analysis, default is True indicating only a json and html report will be generated
            dedup: bool
                a designation of wheather or not to enable deduplication to drop the duplicated reads/pairs, default is False
            threads: int
                an integer representing the number of threads utilized for the operation, default is 1
        """
        self.read_set = read_set
        self.out_dir = out_dir
        self.out_prefix = out_prefix
        self.out_prefix_2 = f"{self.out_prefix}_2"
        self.qualified_quality_phred = qualified_quality_phred
        self.min_read_len = min_read_len
        self.max_read_len = max_read_len
        self.trim_front_bp = trim_front_bp
        self.trim_tail_bp = trim_tail_bp
        self.report_only = report_only
        self.dedup = dedup
        self.threads = threads
        self.paired = self.read_set.is_paired
        

    def run_fastp(self):
        """
        Run the fastp program with the designated paramters selected during the initalization of the class.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        json = os.path.join(self.out_dir,f"{self.out_prefix}.json")
        html = os.path.join(self.out_dir,f"{self.out_prefix}.html")
        out1 = os.path.join(self.out_dir,f"{self.out_prefix}.fastp.fastq")

        self.result_files["html"] = html
        self.result_files["json"] = json

        cmd_args = {'-j':json, '-h':html, '-w':self.threads}
        cmd_args['-i'] = self.read_set.files[0]
        cmd_args['-f'] = self.trim_front_bp
        cmd_args['-t'] = self.trim_tail_bp
        cmd_args['-l'] = self.min_read_len
        cmd_args['--length_limit'] = self.max_read_len
        cmd_args['-q'] = self.qualified_quality_phred
        if self.paired:
            cmd_args['-I'] = self.read_set.files[1]
            out2 = os.path.join(self.out_dir,f"{self.out_prefix_2}.fastp.fastq")
            if not self.report_only:
                cmd_args['-O'] = out2
        if self.dedup:
            cmd_args['-D'] = ''
        if not self.report_only:
            cmd_args['-o'] = out1
            self.result_files["output_files_fastp"].append(out1)
            if self.paired:
                self.result_files["output_files_fastp"].append(out2)
            

        cmd = "fastp {}".format((" ".join(f'{k} {v}' for k,v in cmd_args.items())))
        (self.stdout, self.stderr) = run_command(cmd)
        self.status = self.check_files([json, html] + self.result_files["output_files_fastp"])
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