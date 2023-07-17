#!/usr/bin/env python
from sequenoscope.utils.__init__ import run_command
from sequenoscope.constant import DefaultValues
import os


class KatRunner:
    input_path = None
    ref_path = None
    out_path = None
    result_files = {"sect":{"cvg":"", "tsv":""}, "filter":{"jf27":"", "filtered_fastq":[]}, "hist":{"png_file":"", "json_file":""}}
    error_messages = None
    status = False
    threads = 1
    kmersize = 27

    def __init__(self, input_path, ref_path, out_path, out_prefix, threads = 1, kmersize = DefaultValues.kat_hist_kmer_size):
        """
        Initalize the class with input path, ref_path, and output path

        Arguments:
            input_path: sequence object
                an object that contains the list of sequence files for analysis
            ref_path: str
                a string to the path of reference sequence file
            out_path: str
                a string to the path where the output files will be stored
            threads: int
                an integer representing the number of threads utilized for the operation, default is 1
            kmersize: int
                an integer representing the kmer size utilized for the kat filter method, default is 27
        """
        self.input_path = input_path
        self.ref_path = ref_path
        self.out_path = out_path
        self.out_prefix = out_prefix
        self.threads = threads
        self.kmersize = kmersize

    def kat_sect(self):
        """
        Run the kat sect command on the input sequences and reference sequence.

        Returns: 
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        input_fastq = self.input_path.out_files
        ref_fasta = self.ref_path
        out_file_sect = os.path.join(self.out_path, f"{self.out_prefix}")
        cvg_file = os.path.join(self.out_path, f"{self.out_prefix}-counts.cvg")
        tsv_file = os.path.join(self.out_path, f"{self.out_prefix}-stats.tsv")

        self.result_files["sect"]["cvg"] = cvg_file
        self.result_files["sect"]["tsv"] = tsv_file

        kat_sect_cmd = f"kat sect -t {self.threads} -o {out_file_sect} {ref_fasta} {input_fastq}"
        (self.stdout, self.stderr) = run_command(kat_sect_cmd)
        self.status = self.check_files([cvg_file, tsv_file])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

    def kat_filter(self, exclude=False):
        """
        Run the kat filter kmer command on the reference sequence to generate a kmer hash table and kat filter seq commmand
        on the input sequences along with the kat filter kmer hash table.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        ref_fasta = self.ref_path
        out_file_hash = os.path.join(self.out_path, f"{self.out_prefix}_hash")
        jf_file = os.path.join(self.out_path, f"{self.out_prefix}_hash-in.jf{self.kmersize}")
        paired = self.input_path.is_paired

        self.result_files["filter"][f"jf{self.kmersize}"] = jf_file

        hash_build_command = f"kat filter kmer -m {self.kmersize} -o {out_file_hash} {ref_fasta}"
        (self.stdout, self.stderr) = run_command(hash_build_command)
        self.status = self.check_files([jf_file])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

        out_file_filter = os.path.join(self.out_path, f"{self.out_prefix}_filtered")

        cmd = [ 'kat','filter','seq','-t',f'{self.threads}','-o',out_file_filter,'--seq',self.input_path.files[0]]

        if paired:
            cmd.append('--seq2')
            cmd.append(self.input_path.files[1])
            self.result_files["filter"]["filtered_fastq"].append(os.path.join(self.out_path, f"{self.out_prefix}_filtered.in.R1.fastq"))
            self.result_files["filter"]["filtered_fastq"].append(os.path.join(self.out_path, f"{self.out_prefix}_filtered.in.R2.fastq"))
        else:
            self.result_files["filter"]["filtered_fastq"].append(os.path.join(self.out_path, f"{self.out_prefix}_filtered.in.fastq"))

        if exclude == True:
            cmd.insert(3, "-i")
        
        cmd.append(jf_file)
        cmd_string = " ".join(cmd)

        (self.stdout, self.stderr) = run_command(cmd_string)
        self.status = self.check_files([jf_file] + self.result_files["filter"]["filtered_fastq"])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty, check error message\n{}".format(self.stderr)
            raise ValueError(str(self.error_messages))

    def kat_hist(self):
        """
        Run the kat hist command on the input file sequences

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        input_fastq = self.input_path.out_files
        out_file_hist = os.path.join(self.out_path, f"{self.out_prefix}_histogram_file")
        png_file = os.path.join(self.out_path, f"{self.out_prefix}_histogram_file.png")
        json_file = os.path.join(self.out_path, f"{self.out_prefix}_histogram_file.dist_analysis.json")

        self.result_files["hist"]["png_file"] = png_file
        self.result_files["hist"]["json_file"] = json_file

        kat_hist_cmd = f"kat hist -t {self.threads} -m {self.kmersize} -o {out_file_hist} {input_fastq}"
        (self.stdout, self.stderr) = run_command(kat_hist_cmd)
        self.status = self.check_files([png_file, json_file])
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
