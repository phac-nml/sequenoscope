#!/usr/bin/env python

import os
from math import log
from sequenoscope.constant import DefaultValues
from sequenoscope.utils.parser import fastq_parser
from sequenoscope.analyze.bam import BamProcessor
from sequenoscope.utils.__init__ import is_non_zero_file


class SeqManifest:
    # Fields for the manifest output
    fields = [
        'sample_id', 'read_id', 'read_len', 'read_qscore', 'channel',
        'start_time', 'end_time', 'decision', 'fastp_status',
        'is_mapped', 'is_uniq', 'contig_id',
    ]

    def __init__(self, sample_id, in_bam, out_prefix, out_dir, min_coverage,
                 in_fastq=None, fastp_fastq=None, in_seq_summary=None, read_list=None,
                 start_time=None, end_time=None, delim="\t"):
        """
        Initialize the SeqManifest object with sample and file information.
        
        Raises:
            ValueError: if required sequencing summary or fastq inputs are missing.
        """
        self.delim = delim
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.sample_id = sample_id
        self.in_seq_summary = in_seq_summary
        self.in_fastq = in_fastq
        self.fastp_fastq = fastp_fastq
        self.read_list = read_list
        self.start_time = start_time
        self.end_time = end_time
        self.min_coverage = min_coverage
        self.filtered_reads = {}
        self.raw_reads = {}
        self.status = False
        self.error_messages = None

        if self.in_seq_summary is None:
            if self.start_time is None or self.end_time is None:
                raise ValueError('No sequencing summary specified; please specify a start and end datetime.')
            if self.in_fastq is None:
                raise ValueError('No sequencing summary specified; please provide the initial fastq file for calculations.')

        self.bam_obj = BamProcessor(input_file=in_bam, min_coverage=self.min_coverage)

        if self.fastp_fastq:
            self.process_fastq(self.fastp_fastq, self.filtered_reads)
        if self.in_fastq:
            self.process_fastq(self.in_fastq, self.raw_reads)

        if self.in_seq_summary:
            if not is_non_zero_file(self.in_seq_summary):
                raise ValueError(f"Specified sequencing summary file {self.in_seq_summary} does not exist or is empty.")
            self.create_manifest_with_sum()
        else:
            self.create_manifest_no_sum()

        self.status = True

    @staticmethod
    def error_prob_list_tab(n):
        """Generate a list of error probabilities for qualities 0..n."""
        return [10 ** (q / -10) for q in range(n + 1)]

    def calc_mean_qscores(self, qual, tab=None):
        """Calculate the mean quality score using error probabilities."""
        if tab is None:
            tab = SeqManifest.error_prob_list_tab(DefaultValues.nanoget_threshold)
        if qual:
            avg_error = sum(tab[q] for q in qual) / len(qual)
            return -10 * log(avg_error, 10)
        return 0

    def convert_qscores(self, qual_string):
        """Convert a Phred quality string into a list of integer scores."""
        return [ord(c) - DefaultValues.phred_33_encoding_value for c in qual_string]

    def process_fastq(self, fastq_file_list, read_dict):
        """Process FASTQ files and store read length and computed quality score in a dictionary."""
        for fastq_file in fastq_file_list:
            fastq_obj = fastq_parser(fastq_file)
            for record in fastq_obj.parse():
                read_id = fastq_obj.read_id_from_record
                seq = record[1]
                seq_len = len(seq)
                qual = self.convert_qscores(record[3])
                qscore = self.calc_mean_qscores(qual)
                read_dict[read_id] = [seq_len, qscore]

    def create_row(self):
        """Create an empty row dictionary with keys from fields."""
        return {field: '' for field in self.fields}

    def create_manifest_with_sum(self):
        """Create the manifest file using a sequencing summary."""
        manifest_file = os.path.join(self.out_dir, f"{self.out_prefix}.txt")
        with open(self.read_list, 'r') as file:
            read_set = {line.strip() for line in file if line.strip() != 'read_id'}

        with open(manifest_file, 'w') as fout, open(self.in_seq_summary, 'r') as fin:
            header = next(fin).strip().split(self.delim)
            fout.write("\t".join(self.fields) + "\n")
            for line in fin:
                row = line.strip().split(self.delim)
                if len(row) < len(header):
                    continue
                row_data = dict(zip(header, row))
                read_id = row_data.get('read_id')
                if read_id not in read_set:
                    continue

                read_len = row_data.get('sequence_length_template', 0)
                read_qual = row_data.get('mean_qscore_template', 0)

                is_uniq = True
                is_mapped = False
                start_time_val = row_data.get('start_time', '')
                duration = row_data.get('duration', '')
                if not start_time_val:
                    start_time_val = self.start_time
                    end_time_val = self.end_time
                else:
                    start_time_val = float(start_time_val)
                    end_time_val = start_time_val + float(duration) if duration else ''

                out_row = self.create_row()
                for field in self.fields:
                    if field in row_data:
                        out_row[field] = row_data[field]

                mapped_contigs = [cid for cid, stats in self.bam_obj.ref_stats.items()
                                  if read_id in stats['reads'] and cid != '*']
                if mapped_contigs:
                    is_mapped = True
                    if len(mapped_contigs) > 1:
                        is_uniq = False

                fastp_status = read_id in self.filtered_reads

                out_row.update({
                    'fastp_status': fastp_status,
                    'sample_id': self.sample_id,
                    'read_id': read_id,
                    'is_mapped': is_mapped,
                    'is_uniq': is_uniq,
                    'read_len': read_len,
                    'read_qscore': read_qual,
                    'start_time': start_time_val,
                    'end_time': end_time_val,
                    'decision': row_data.get('end_reason', '')
                })

                if not mapped_contigs:
                    out_row['contig_id'] = ''
                    fout.write("\t".join(str(x) for x in out_row.values()) + "\n")
                else:
                    for contig_id in mapped_contigs:
                        out_row['contig_id'] = contig_id
                        fout.write("\t".join(str(x) for x in out_row.values()) + "\n")

        if not self.check_files([manifest_file]):
            raise ValueError("One or more files were not created or were empty")

    def create_manifest_no_sum(self):
        """Create the manifest file when no sequencing summary is provided, using a read list and raw FASTQ data."""
        manifest_file = os.path.join(self.out_dir, f"{self.out_prefix}.txt")
        with open(manifest_file, 'w') as fout, open(self.read_list, 'r') as fin:
            fout.write("\t".join(self.fields) + "\n")
            header = next(fin).strip().split(self.delim)
            for line in fin:
                row = line.strip().split(self.delim)
                row_data = dict(zip(header, row))
                read_id = row_data.get('read_id')

                read_len, read_qual = (0, 0)
                if read_id in self.raw_reads:
                    read_len, read_qual = self.raw_reads[read_id]

                is_uniq = True
                is_mapped = False
                start_time_val = self.start_time
                end_time_val = self.end_time

                out_row = self.create_row()
                for field in self.fields:
                    if field in row_data:
                        out_row[field] = row_data[field]

                mapped_contigs = []
                for contig_id, stats in self.bam_obj.ref_stats.items():
                    if read_id in stats['reads']:
                        read_len, read_qual = stats['reads'][read_id]
                        if contig_id != '*':
                            mapped_contigs.append(contig_id)
                if mapped_contigs:
                    is_mapped = True
                if len(mapped_contigs) > 1:
                    is_uniq = False

                fastp_status = (read_id in self.filtered_reads)
                out_row.update({
                    'fastp_status': fastp_status,
                    'sample_id': self.sample_id,
                    'read_id': read_id,
                    'is_mapped': is_mapped,
                    'is_uniq': is_uniq,
                    'read_len': read_len,
                    'read_qscore': read_qual,
                    'start_time': start_time_val,
                    'end_time': end_time_val,
                    'decision': "signal_positive",
                    'channel': "1"
                })

                if not mapped_contigs:
                    out_row['contig_id'] = ''
                    fout.write("\t".join(str(x) for x in out_row.values()) + "\n")
                else:
                    for contig_id in mapped_contigs:
                        out_row['contig_id'] = contig_id
                        fout.write("\t".join(str(x) for x in out_row.values()) + "\n")

        if not self.check_files([manifest_file]):
            raise ValueError("One or more files were not created or were empty")

    def check_files(self, files_to_check):
        """Check if each file in the list exists and is non-empty."""
        if isinstance(files_to_check, str):
            files_to_check = [files_to_check]
        for f in files_to_check:
            if not os.path.isfile(f) or os.path.getsize(f) == 0:
                return False
        return True


class SeqManifestSummary:
    # Define the default fields; updated with new dynamic field naming.
    fields = [
        'sample_id', 'est_genome_size', 'est_coverage', 'total_bases', 'total_fastp_bases',
        'mean_read_length', 'taxon_id', 'taxon_length', 'taxon_mean_coverage',
        # Dynamic field for taxon covered bases using min_cov:
        'taxon_covered_bases',  
        'taxon_%_covered_bases', 'total_taxon_ref_mapped_bases', 'taxon_mean_read_length'
    ]

    def __init__(self, sample_id, bam_obj, out_prefix, out_dir, genome_size, coverage,
                 fastp_json_file=None, paired=False):
        self.sample_id = sample_id
        self.fastp_json_file = fastp_json_file
        self.bam_obj = bam_obj
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.genome_size = genome_size
        self.coverage = coverage
        self.paired = paired

        # Dynamic field: taxon_covered_bases_<min_cov>X
        min_cov = self.bam_obj.min_coverage
        self.taxon_coverage_field = f"taxon_covered_bases_{min_cov}X"
        # New dynamic field for taxon percentage covered:
        self.taxon_percentage_field = f"taxon_%_covered_bases_{min_cov}X"
        # New dynamic field for total mapped bases:
        self.total_taxon_ref_mapped_field = "total_taxon_ref_mapped_bases"  # remains fixed

        # Update fields list to include the dynamic names.
        self.fields = [
            'sample_id', 'est_genome_size', 'est_coverage', 'total_bases', 'total_fastp_bases',
            'mean_read_length', 'taxon_id', 'taxon_length', 'taxon_mean_coverage',
            self.taxon_coverage_field,
            self.taxon_percentage_field, self.total_taxon_ref_mapped_field, 'taxon_mean_read_length'
        ]

    def create_row(self):
        """Create an empty row dictionary for the summary manifest."""
        return {field: '' for field in self.fields}

    def generate_summary(self):
        """
        Generate the summary manifest file using BAM and fastp JSON data.
        """
        summary_manifest_file = os.path.join(self.out_dir, f"{self.out_prefix}.txt")
        with open(summary_manifest_file, 'w') as fout:
            fout.write("\t".join(self.fields) + "\n")
            for contig_id, stats in self.bam_obj.ref_stats.items():
                out_row = self.create_row()
                out_row["sample_id"] = self.sample_id
                out_row["est_genome_size"] = self.genome_size
                out_row["est_coverage"] = self.coverage
                out_row["total_bases"] = self.fastp_json_file["summary"]["before_filtering"]["total_bases"]
                out_row["total_fastp_bases"] = self.fastp_json_file["summary"]["after_filtering"]["total_bases"]
                out_row["mean_read_length"] = self.fastp_json_file["summary"]["after_filtering"]["read1_mean_length"]
                out_row["taxon_id"] = contig_id
                out_row["taxon_length"] = stats['length']
                out_row[self.taxon_coverage_field] = stats['covered_bases']
                if stats['length'] != 0:
                    out_row[self.taxon_percentage_field] = (stats['covered_bases'] / stats['length']) * 100
                    out_row[self.total_taxon_ref_mapped_field] = stats['total_mapped_bases']
                    out_row["taxon_mean_read_length"] = stats['mean_len']
                    out_row["taxon_mean_coverage"] = stats['mean_cov']
                else:
                    out_row[self.taxon_percentage_field] = 0
                    out_row[self.total_taxon_ref_mapped_field] = 0
                    out_row["taxon_mean_read_length"] = 0
                    out_row["taxon_mean_coverage"] = 0
                fout.write("\t".join(str(x) for x in out_row.values()) + "\n")

        if not self.check_files([summary_manifest_file]):
            raise ValueError("One or more files were not created or were empty")

    def check_files(self, files_to_check):
        """Check if the given file(s) exist and are non-empty."""
        if isinstance(files_to_check, str):
            files_to_check = [files_to_check]
        for f in files_to_check:
            if not os.path.isfile(f) or os.path.getsize(f) < 0:
                return False
        return True