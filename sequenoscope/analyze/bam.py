#!/usr/bin/env python

import statistics
import pysam
from math import log
from sequenoscope.constant import DefaultValues
from sequenoscope.utils.__init__ import run_command, is_non_zero_file



class BamProcessor:

    alignment_file = None
    index_file = None
    pysam_obj = None
    threads = 1
    read_locations = {}
    ref_stats = {}
    ref_coverage = {}
    status = True
    error_msg = ''

    def __init__(self,input_file, min_coverage):
        """
        Initalize the class with an input bam file

        Arguments:
            input_file: str
                a string that designates the path of the bam file to be analyzed
        """
        self.alignment_file = input_file
        self.min_coverage = min_coverage
        if not is_non_zero_file(input_file):
            self.status = False
            self.error_msg = "Error bam file {} does not exist".format(input_file)
            return
        index_file = "{}.bai".format(input_file)
        self.index_file = index_file
        if not is_non_zero_file(index_file):
            (stdout,stderr) =self.index_bam()
        if not is_non_zero_file(index_file):
            self.status = False
            self.error_msg = "STDOUT:{}\nSTDERR:{}".format(stdout,stderr)
            return
        self.ref_stats = self.get_bam_stats()
        self.init_base_cov()
        self.pysam_obj = pysam.AlignmentFile(input_file, "rb")
        self.process_bam()


    def init_base_cov(self):
        """
        Uses the contig lengths from self.ref_stats to create a list of positions with counts initialized to 0
        """
        for contig_id in self.ref_stats:
            self.ref_coverage[contig_id] = [0] * self.ref_stats[contig_id]['length']


    def process_bam(self):
        """
        Reads a bam file line by line and produces summary statistics based on each contig
        """
        for contig_id in self.ref_stats:
            contig_len = self.ref_stats[contig_id]['length']
            #covered_positions = set()
            num_reads = 0
            total_bases = 0
            lengths = []
            qualities = []
            for read in self.pysam_obj.fetch(contig_id):
                num_reads+=1
                read_id = read.query_name
                seq = read.query_sequence
                if seq is not None:
                    length = len(seq)
                else:
                    length = 0
                total_bases += length
                lengths.append(length)
                qual = read.query_qualities
                qscore = self.calc_mean_qscores(qual)
                qualities.append(qscore)
                self.ref_stats[contig_id]['reads'][read_id] = (length,qscore)
                if contig_id == '*':
                    continue
                start_pos = read.reference_start
                aln_len = read.query_alignment_length
                for i in range(start_pos,start_pos+aln_len):
                    if i < contig_len: #and i not in covered_positions:
                        self.ref_coverage[contig_id][i]+=1
                        #covered_positions.add(i)

            lengths = sorted(lengths,reverse=True)
            if len(self.ref_coverage[contig_id]) > 0:
                self.ref_stats[contig_id]['mean_cov'] = statistics.mean(self.ref_coverage[contig_id])
                self.ref_stats[contig_id]['covered_bases'] = self.count_cov_bases(self.ref_coverage[contig_id])
                self.ref_stats[contig_id]['total_mapped_bases'] = sum(self.ref_coverage[contig_id])
            self.ref_stats[contig_id]['n50'] = self.calc_n50(lengths,total_bases)
            self.ref_stats[contig_id]['num_reads'] = num_reads
            if len(lengths) > 0:
                self.ref_stats[contig_id]['median_len'] = statistics.median(lengths)
                self.ref_stats[contig_id]['mean_len'] = statistics.mean(lengths)
            if len(qualities) > 0:
                self.ref_stats[contig_id]['median_qual'] = statistics.median(qualities)
                self.ref_stats[contig_id]['mean_qual'] = statistics.mean(qualities)
        return

    def calc_n50(self,lengths,total_length):
        """
        Calculates the N50 of a set of read lengths

        Arguments:
            lengths: list
                a list of read lengths for N50 calcualtion
            total_length: int
                total number of bases accross all reads
        
        Returns: 
            int:
               tabulated N50 value based on lengths 
        """
        target_len = int(total_length / 2)
        s = 0
        try:
            global l
            for l in lengths:
                if s >= target_len:
                    return l
                s+=l
            return l
        except NameError:
            l = 0
            return l

    def count_cov_bases(self,list_of_values, min_value=None, max_value=9999999999999):
        """
        Counts positions where the count is >=min and <= max

        Arguments:
            list_of_values: list
                list of coverage value for calcualtion
            min_value: int
                minimum coverage value
            max_value: int
                maximum coverage value

        Returns:
            int:
                int number of positions meeting this threshold
        """
        if min_value is None:
            min_value = self.min_coverage
        
        total = 0
        for v in list_of_values:
            if v >= min_value and v<=max_value:
                total+=1
        return total
    
    def error_prob_list_tab(n):
        """
        generate a list of error rates for qualities less than or equal to n.
        source: github.com/wdecoster/nanoget/blob/master/nanoget/utils.py

        Arguments: 
            n: error probability threshold

        Returns:
            list:
                list of error rates
        """
        return [10**(q / -10) for q in range(n+1)]

    def calc_mean_qscores(self,qual,tab=error_prob_list_tab(DefaultValues.nanoget_threshold)):
        """
        Calculates the mean quality score for a read where they have been converted to Phred.
        Phred scores are first converted to probabilites, then the average error probability is calculated.
        The average is then converted back to the Phred scale.

        Arguments:
            qual: string
                string of Phred 33 ints for quality calcualtions
            
            tab: list
                list of error rates for qaulties specified
                default designation of 128 based on nanoget 

        Returns:
            float:
                mean qscore
        """
        if qual:
            phred_score = -10 * log(sum([tab[q] for q in qual]) / len(qual) , 10)
            return phred_score
        else:
            return 0


    def get_bam_stats(self):
        """
        Wrapper class around SAMTOOLS IDXSTATS for getting information about bam file

        Returns:
            dictionary:
                dictionary with some results from the samtools IDXSTATS tool

        """
        cmd = [
            'samtools',
            'idxstats',
            "{}".format(self.alignment_file)
        ]
        cmd = " ".join(cmd)
        (stdout,stderr) = run_command(cmd)
        result = {}
        stdout = stdout.split("\n")
        for row in stdout:
            row = row.split("\t")
            if len(row) < DefaultValues.samtools_idxstats_field_number:
                continue
            result[row[0]] = {'length':int(row[1]),
                              'reads': {},'num_reads':0,'mean_cov':0,
                              'covered_bases':0,'mean_len':0,'median_len':0,
                              'mean_qual':0,'median_qual':0,'n50':0}
        return result


    def index_bam(self):
        """
        Wrapper class around SAMTOOLS INDEX for creating a bam index file

        Returns:
            file object:
                indexed bam file
        """
        cmd = [
            'samtools',
            'index',
            '-@',"{}".format(self.threads),
            "{}".format(self.alignment_file),
            "{}".format(self.index_file)
        ]
        cmd = " ".join(cmd)
        return run_command(cmd)