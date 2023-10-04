#!/usr/bin/env python
import os
from sequenoscope.constant import DefaultValues

class FastqExtractor:
    out_prefix = None
    out_dir = None
    read_set = None
    status = False
    result_files = {"read_list_file":""}
    
    def __init__(self, read_set, out_prefix, out_dir):
        """
        Initalize the class with read_set, out_prefix, and out_dir

        Arguments:
            read_set: sequence object
                an object that contains the list of sequence files for analysis
            out_prefix: str
                a designation of what the output files will be named
            out_dir: str
                a string to the path where the output files will be stored
        """
        self.reads = []
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.read_set = read_set
   
    def extract_single_reads(self):
        """
        Extracts the read ids from a single-end fastq file based on the paramters intialized in the previous method.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        self.extract_reads(self.read_set.files[0], self.reads, read_delimiter=" ", split_delimiter=None)

        self.write_reads(read_lists=[self.reads])
               
    def extract_paired_reads(self):
        """
        Extracts the read ids from a paired-end fastq file based on the paramters intialized in the previous method.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        forward_reads = []
        self.extract_reads(self.read_set.files[0], forward_reads, split_delimiter=":")
        reverse_reads = []
        self.extract_reads(self.read_set.files[1], reverse_reads, split_delimiter=":")

        self.write_reads(read_lists=[forward_reads, reverse_reads])
        
    def alt_extract_paired_reads(self):
        """
        Extracts the read ids from a paired-end fastq file based on the paramters intialized in the previous method.

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        forward_reads = []
        with open(self.read_set.files[0], 'r') as f:
            for line in f:
                if line.startswith('@'):
                    if line.endswith('1\n'):
                        read_id = line.strip().split()[0][1:] #+ '_R1'
                        forward_reads.append(read_id)
        reverse_reads = []
        with open(self.read_set.files[1], 'r') as f:
            for line in f:
                if line.startswith('@'):
                    if line.endswith('2\n'):
                        read_id = line.strip().split()[0][1:] #+ '_R2'
                        reverse_reads.append(read_id)
        self.write_reads(read_lists=[forward_reads, reverse_reads])
    
    def extract_reads(self, file, read_list, read_delimiter=None, split_delimiter=None):
        """
        Strip the lines of fastq file based on various delimitors that result
        from different sequencing instrumnetaion output

        Arguments:
            file: list object
                file in read set for extracting reads
            read_list: list
                list of where to store reads
            read_delimitor: str
                delimitor that is located by the read id after stripping the lines
            split_delimitor: str
                delimitor that is located by the read id before stripping the lines
        """
        line_count = 0
        with open(file, 'r') as f:
            for line in f:
                line_count += 1
                if line.startswith(DefaultValues.fastq_line_starter) and (line_count % 4 == 1):
                    if len(line.strip().split(split_delimiter)) >= DefaultValues.fastq_sample_row_number:
                        read_id = line.strip().split(read_delimiter)[0][1:]
                        read_list.append(read_id)
                    else:
                        read_id = line.strip().split()[0][1:]
                        read_list.append(read_id)
    
    def write_reads(self, read_lists=[]):
        """
        Append the lines extracted from a list into a file output and check if the file was created

        Arguments:
            read_lists: list
                lists of where reads are stored

        Returns:
            bool:
                returns True if the generated output file is found and not empty, False otherwise
        """
        output_file = os.path.join(self.out_dir,f"{self.out_prefix}.txt")
        self.result_files["read_list_file"] = output_file

        if len(read_lists) == 1:
            with open(output_file, 'w') as f:
                f.write("read_id\n")  # Write the header row
                for read in read_lists[0]:
                    f.write(f"{read}\n")
        else:
            with open(output_file, 'w') as f:
                f.write("read_id\n")  # Write the header row
                for read in read_lists[0] + read_lists[1]:
                    f.write(f"{read}\n")
            
        self.status = self.check_files(output_file)
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