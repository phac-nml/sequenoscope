#!/usr/bin/env python

from dataclasses import dataclass

@dataclass(frozen=True)
class SequenceTypes:
    paired_end: str = 'PE'
    single_end: str = 'SE'

@dataclass(frozen=True)
class DefaultValues:
    minimap2_kmer_size: int = 15
    kat_hist_kmer_size: int = 27
    nanoget_threshold: int = 128
    samtools_idxstats_field_number: int = 4
    fastq_sample_row_number: int = 4
    fastq_line_starter: str = "@"
    phred_33_encoding_value: int = 33
    max_nanopore_channel: int = 512