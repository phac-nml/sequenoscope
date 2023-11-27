#!/usr/bin/env python

import os
from math import log
from sequenoscope.constant import DefaultValues
from sequenoscope.utils.parser import fastq_parser
from sequenoscope.analyze.bam import BamProcessor
from sequenoscope.utils.__init__ import is_non_zero_file


class SeqManifest:
    fields = [
        'sample_id','read_id','read_len','read_qscore','channel',
        'start_time','end_time','decision','fastp_status',
        'is_mapped','is_uniq','contig_id',
    ]
    sample_id = ''
    in_seq_summary = ''
    read_list = ''
    out_prefix = ''
    out_dir = ''
    bam_obj = None
    fastp_obj = None
    fastp_fastq = []
    delim = "\t"
    status = False
    start_time = ''
    end_time = ''
    filtered_reads = {}
    raw_reads = {}
    error_messages = None

    def __init__(self,sample_id,in_bam,out_prefix, out_dir, min_coverage, in_fastq=None, fastp_fastq=None,in_seq_summary=None,
                  read_list=None,start_time=None,end_time=None,delim="\t"):
        """
        Initalize the class with sample_id, in_bam, out_prefix, and out_dir. Analyze reads based on seq summary and 
        fastp fast availbility by producing manifest files.

        Arguments:
            sample_id: str
                a string of the name of the sample to be analyzed
            in_bam: str
                a string to the path where the bam file is stored
            out_prefix: str
                a designation of what the output files will be named
            out_dir: str
                a designation of the output directory
            fastq_fastq: str
                a designation of where the filitered fastq produced by fastp is stored.
            in_seq_summary: str
                a designation of where the sequencing summary produced by the Nanopore sequencers is stored
            read_list: str
                a designation of where the read list produced from the original fastq is stored
            start_time: int
                an integer representing the start time when seq summary isn't provided.
            end_time: int
                an integer representing the end time when seq summary isn't provided.
            delim: str
                a string that designates the delimiter used to parse files. default is tab delimiter
        """
        self.delim = delim
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.in_fastq = in_fastq
        self.sample_id = sample_id
        self.in_seq_summary = in_seq_summary
        self.fastp_fastq = fastp_fastq
        self.start_time = start_time
        self.end_time = end_time
        self.read_list = read_list
        self.min_coverage = min_coverage

        if self.in_seq_summary is None:
            if self.start_time is None or self.end_time is None:
                self.status = False
                self.error_msg = 'Error no sequence summary specified, please specify a start and end datetime'
                return
            if self.in_fastq is None:
                self.status = False
                self.error_msg = 'Error no sequence summary specified, please add a the intial fastq file for calculations'
                return

        self.bam_obj = BamProcessor(input_file=in_bam, min_coverage=self.min_coverage)

        if self.fastp_fastq is not None:
            self.process_fastq(self.fastp_fastq, self.filtered_reads)
        if self.in_fastq is not None: 
            self.process_fastq(self.in_fastq, self.raw_reads)

        if self.in_seq_summary is not None and not is_non_zero_file(self.in_seq_summary):
            self.status = False
            self.error_msg = f"Error specified seq summary file {self.in_seq_summary} does not exist"
            return
        if self.in_seq_summary is not None:
            self.create_manifest_with_sum()
        else:
            self.create_manifest_no_sum()
    
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

        Returns:
            float:
                mean qscore
        """
        if qual:
            phred_score = -10 * log(sum([tab[q] for q in qual]) / len(qual) , 10)
            return phred_score
        else:
            return 0


    def convert_qscores(self,qual_string):
        """
        Calculates the mean quality score for a read where they have been converted to Phred

        Arguments
            qual_string: string of phred 33 ints for quality

        Returns:
            float:
                mean qscore
        """
        qual_values = []
        for c in qual_string:
            qual_values.append(ord(c) - DefaultValues.phred_33_encoding_value)
        return qual_values

    def process_fastq(self, fastq_file_list, read_dict):
        """
        Process the fastq file and extract reads, quality, and qscores

        Argument:
            fastq_file_list:
                list of fastq files
            read_dict:
                dictonary to store reads
        """
        for fastq_file in fastq_file_list:
            fastq_obj = fastq_parser(fastq_file)
            for record in fastq_obj.parse():
                read_id = fastq_obj.read_id_from_record
                seq = record[1]
                seq_len = len(seq)
                qual = self.convert_qscores(record[3])
                qscore = self.calc_mean_qscores(qual)
                read_dict[read_id] = [seq_len,qscore]


    def create_row(self):
        """
        create rows and store them into a dictionary

        Returns:
            dict:
                dictionary of rows produced.
        """
        out_row = {}
        for field_id in self.fields:
            out_row[field_id] = ''
        return out_row
        


    def create_manifest_with_sum(self):
        """
        Create a manifest file with various statistics when a sequencing summary is present

        Returns: 
            bool: 
                True if the summary manifest file was created, False otherwise.
        """
        manifest_file = os.path.join(self.out_dir,f"{self.out_prefix}.txt")
        fout = open(manifest_file,'w')
        fout.write("{}\n".format("\t".join(self.fields)))

        fin = open(self.in_seq_summary,'r')
        header = next(fin).strip().split(self.delim)

        read_set = set()
        with open(self.read_list, 'r') as file:
            lines = file.readlines()

        for line in lines:
            line = line.strip()
            if line != 'read_id':
                read_set.add(line)

        for line in fin:
            row = line.strip().split(self.delim)
            row_data = {}
            for i in range(0,len(row)) :
                row_data[header[i]] = row[i]

            read_id = row_data['read_id']

            if read_id not in read_set:
              continue
            try:
                read_len = row_data['sequence_length_template']
                read_qual = row_data['mean_qscore_template']
            except KeyError:
                read_len = 0
                read_qual = 0
            
            is_uniq = True
            is_mapped = False
            start_time = row_data['start_time']
            end_time = ''
            duration = row_data['duration']
            if start_time == '':
                start_time = self.start_time
                end_time = self.end_time
            else:
                start_time = float(start_time)
                if duration != '':
                    end_time = start_time + float(duration)

            out_row = self.create_row()
            for field_id in self.fields:
                if field_id in row_data:
                    out_row[field_id] = row_data[field_id]
            mapped_contigs = []
            for contig_id in self.bam_obj.ref_stats:
                if read_id in self.bam_obj.ref_stats[contig_id]['reads']:
                    pass
                    if contig_id != '*':
                        mapped_contigs.append(contig_id)

            if len(mapped_contigs) > 0:
                is_mapped = True
            if len(mapped_contigs) > 1:
                is_uniq = False
            

            fastp_status = False
            if read_id in self.filtered_reads:
                fastp_status = True

            out_row['fastp_status'] = fastp_status
            out_row['sample_id'] = self.sample_id
            out_row['read_id'] = read_id
            out_row['is_mapped'] = is_mapped
            out_row['is_uniq'] = is_uniq
            out_row['read_len'] = read_len
            out_row['read_qscore'] = read_qual
            out_row['start_time'] = start_time
            out_row['end_time'] = end_time
            out_row['decision'] = row_data['end_reason']

            if len(mapped_contigs) == 0:
                out_row['contig_id'] = ''
                fout.write("{}\n".format("\t".join([str(x) for x in out_row.values()])))

            for contig_id in mapped_contigs:
                out_row['contig_id'] = contig_id
                fout.write("{}\n".format("\t".join([str(x) for x in out_row.values()])))

        self.status = self.check_files([manifest_file])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty"
            raise ValueError(str(self.error_messages))

        fin.close()
        fout.close()

    def create_manifest_no_sum(self):
        """
        Create a manifest file with various statistics when a sequencing summary is NOT present. Uses read list 
        file instead.

        Returns: 
            file object: 
                seq manifest text file
        """
        
        manifest_file = os.path.join(self.out_dir,f"{self.out_prefix}.txt")
        fout = open(manifest_file,'w')
        fout.write("{}\n".format("\t".join(self.fields)))

        fin = open(self.read_list,'r')
        header = next(fin).strip().split(self.delim)

        for line in fin:
            row = line.strip().split(self.delim)
            row_data = {}
            for i in range(0,len(row)) :
                row_data[header[i]] = row[i]

            read_id = row_data['read_id']

            read_len = 0
            read_qual = 0
            if read_id in self.raw_reads:
                read_len = self.raw_reads[read_id][0]
                read_qual = self.raw_reads[read_id][1]

            is_uniq = True
            is_mapped = False
            start_time = self.start_time
            end_time = self.end_time

            out_row = self.create_row()
            for field_id in self.fields:
                if field_id in row_data:
                    out_row[field_id] = row_data[field_id]
            mapped_contigs = []
            for contig_id in self.bam_obj.ref_stats:
                if read_id in self.bam_obj.ref_stats[contig_id]['reads']:
                    read_len = self.bam_obj.ref_stats[contig_id]['reads'][read_id][0]
                    read_qual = self.bam_obj.ref_stats[contig_id]['reads'][read_id][1]
                    if contig_id != '*':
                        mapped_contigs.append(contig_id)

            if len(mapped_contigs) > 0:
                is_mapped = True
            if len(mapped_contigs) > 1:
                is_uniq = False

            fastp_status = False
            if read_id in self.filtered_reads:
                fastp_status = True
            out_row['fastp_status'] = fastp_status
            out_row['sample_id'] = self.sample_id
            out_row['read_id'] = read_id
            out_row['is_mapped'] = is_mapped
            out_row['is_uniq'] = is_uniq
            out_row['read_len'] = read_len
            out_row['read_qscore'] = read_qual
            out_row['start_time'] = start_time
            out_row['end_time'] = end_time
            out_row['decision'] = "signal_positive"
            out_row['channel'] = "1"


            if len(mapped_contigs) == 0:
                out_row['contig_id'] = ''
                fout.write("{}\n".format("\t".join([str(x) for x in out_row.values()])))

            for contig_id in mapped_contigs:
                out_row['contig_id'] = contig_id
                fout.write("{}\n".format("\t".join([str(x) for x in out_row.values()])))

        self.status = self.check_files([manifest_file])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty"
            raise ValueError(str(self.error_messages))
        
        fin.close()
        fout.close()

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
    
class SeqManifestSummary:
    fields = [
        'sample_id','est_genome_size','est_coverage','total_bases','total_fastp_bases',
        'mean_read_length','taxon_id','taxon_length','taxon_covered_bases',
        'taxon_%_covered_bases','taxon_mean_read_length'
    ]
    sample_id = ''
    out_prefix = ''
    out_dir = ''
    genome_size = None
    coverage = None
    bam_obj = None
    status = False
    error_messages = None

    def __init__(self,sample_id, bam_obj, out_prefix, out_dir, genome_size, coverage,
                  fastp_json_file=None, paired=False):
        """
        Initalize the class with sample_id, bam_obj, out_prefix, and out_dir. Extract sequencing
        statisitics from various files and append to a summary file

        Arguments:
            sample_id: str
                a string of the name of the sample to be analyzed
            bam_obj: str
                an object that of the SeqManifest class that stores bam information
            out_prefix: str
                a designation of what the output files will be named
            out_dir: str
                a designation of the output directory
            genome_size: int
                a designation of genome size based on mash distance calcualtion.
            coverage: int
                a designation of coverage based on mash distance calcualtion.
            fastp_json_file: str
                a designation to the path of the json file generated from fastp.
            paired: bool
                a designation of wheather or not the files specified belong to paired-end sequencing data
        """
        self.sample_id = sample_id
        self.fastp_json_file = fastp_json_file
        self.bam_obj = bam_obj
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.genome_size = genome_size
        self.coverage = coverage
        self.paired = paired

        # Define the dynamic field name and store it as an instance attribute
        min_cov = self.bam_obj.min_coverage
        self.taxon_coverage_field = f"taxon_covered_bases_{min_cov}X"

        # Use the dynamic field name in the fields list
        self.fields = [
            'sample_id', 'est_genome_size', 'est_coverage', 'total_bases', 'total_fastp_bases',
            'mean_read_length', 'taxon_id', 'taxon_length', 'taxon_mean_coverage',
            self.taxon_coverage_field,  # Using the dynamic field name
            'taxon_%_covered_bases', 'total_taxon_mapped_bases', 'taxon_mean_read_length' 
        ]
    
    def create_row(self):
        """
        create rows and store them into a dictionary

        Returns:
            dict:
                dictionary of rows produced.
        """
        out_row = {}
        for field_id in self.fields:
            out_row[field_id] = ''
        return out_row
    
    def generate_summary(self):
        """
        Create a summary manifest file with various statistics from different file sources.

        Returns: 
            bool: 
                True if the summary manifest file was created, False otherwise.
        """
        summary_manifest_file = os.path.join(self.out_dir,f"{self.out_prefix}.txt")
        fout = open(summary_manifest_file,'w')
        fout.write("{}\n".format("\t".join(self.fields)))

        out_row = self.create_row()
        for contig_id in self.bam_obj.ref_stats:
            out_row["sample_id"] = self.sample_id
            out_row["est_genome_size"] = self.genome_size
            out_row["est_coverage"] = self.coverage
            out_row["total_bases"] = self.fastp_json_file["summary"]["before_filtering"]["total_bases"]
            out_row["total_fastp_bases"] = self.fastp_json_file["summary"]["after_filtering"]["total_bases"]
            out_row["mean_read_length"] = self.fastp_json_file["summary"]["after_filtering"]["read1_mean_length"]
            out_row["taxon_id"] = contig_id
            out_row["taxon_length"] = self.bam_obj.ref_stats[contig_id]['length']
            out_row[self.taxon_coverage_field] = self.bam_obj.ref_stats[contig_id]['covered_bases']
            if self.bam_obj.ref_stats[contig_id]['length'] != 0:
                out_row["taxon_%_covered_bases"] = ((self.bam_obj.ref_stats[contig_id]['covered_bases']/self.bam_obj.ref_stats[contig_id]['length']) * 100)
                out_row["total_taxon_mapped_bases"] = self.bam_obj.ref_stats[contig_id]['total_mapped_bases']
                out_row["taxon_mean_read_length"] = self.bam_obj.ref_stats[contig_id]['mean_len']
                out_row["taxon_mean_coverage"] = self.bam_obj.ref_stats[contig_id]['mean_cov']
            else:
                out_row["taxon_%_covered_bases"] = 0
                out_row["total_taxon_mapped_bases"] = 0
                out_row["taxon_mean_read_length"] = 0
                out_row["taxon_mean_coverage"] = 0
                
            fout.write("{}\n".format("\t".join([str(x) for x in out_row.values()])))
        
        self.status = self.check_files([summary_manifest_file])
        if self.status == False:
            self.error_messages = "one or more files was not created or was empty"
            raise ValueError(str(self.error_messages))
        
        fout.close()

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
            elif os.path.getsize(f) < 0:
                return False
        return True