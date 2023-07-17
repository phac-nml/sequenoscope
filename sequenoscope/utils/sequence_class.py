#!/usr/bin/env python

class Sequence:
    technology = None
    files = []
    is_paired = False
    out_files = ''

    def __init__(self, tech_name, list_of_seq):
        self.technology = tech_name
        self.files = list_of_seq
        for file in self.files:
            if not self.is_fastq(file):
                raise ValueError(f'{file} is not a valid fastq file')
        self.classify_seq()
        self.output_formatted_files()
        return

    def classify_seq(self):
        if len(self.files) == 2:
            self.is_paired = True

    def is_fastq(self, input):
        with open(input, "r") as f:
            first_line = f.readline().strip()
            if not first_line.startswith("@"):
                return False

            second_line =  f.readline().strip()

            third_line = f.readline().strip()
            if not third_line.startswith("+"):
                return False
            
            fourth_line = f.readline().strip()

            if len(second_line) != len(fourth_line):
                return False
            
        return True

    def output_formatted_files(self):
        self.out_files = ' '.join(self.files)