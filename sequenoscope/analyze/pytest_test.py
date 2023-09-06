#!/usr/bin/env python
from sequenoscope.utils.sequence_class import Sequence
from sequenoscope.utils.parser import GeneralSeqParser
from sequenoscope.analyze.kat import KatRunner
from sequenoscope.analyze.fastP import FastPRunner
from sequenoscope.analyze.minimap2 import Minimap2Runner
from sequenoscope.analyze.processing import SamBamProcessor
from sequenoscope.analyze.bam import BamProcessor
from sequenoscope.analyze.seq_manifest import SeqManifest
from sequenoscope.analyze.fastq_extractor import FastqExtractor
from sequenoscope.utils.parser import FastqPairedEndRenamer
from sequenoscope.analyze.seq_manifest import SeqManifestSummary
from sequenoscope.analyze.mash import MashSketcher

path_ref_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/lambda_genome_reference.fasta"
path_enriched_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Test_br1_sal_lam_enriched.fastq"
path_control_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Test_br1_sal_lam_control.fastq"
path_output = "/home/ameknas/sequenoscope-1/sequenoscope/analyze"
technology = "ONT"

invalid_ref_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/invalid_reference.fasta"
invalid_seq_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/invalid_input.fastq"

tsv_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/output-stats.tsv"
json_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/histogram_file.dist_analysis.json"
sam_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/test_output.sam"
bam_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/test_output.bam"
AS_report_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/adaptive_sampling_report.csv"
UB_csv_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/test_output_unblocked_ids.csv"
SR_csv_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/test_output_stop_receiving.csv"
ND_csv_test_file = "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/test_output_no_decision.csv"
cmd = "python -m sequenoscope.main analzye --input_fastq /home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Test_br1_sal_lam_enriched.fastq --input_reference /home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/invalid_reference.fasta -o test -seq_type sr -min_len 3 -trm_tail 100"
test_bam = "/home/ameknas/sequenoscope-1/test/sample_mapped_bam.bam"

def test_make_test():
    print ("hello world")
    pass

# def test_kat_sect():
#     enriched_sample = Sequence(technology, [path_enriched_test_file, path_enriched_test_file])
#     kat_run = KatRunner(enriched_sample, path_ref_file, path_output, "test")
#     kat_run.kat_sect()
#     assert kat_run.status == True
#     pass

# def test_kat_filter():
#     enriched_sample = Sequence(technology, [path_enriched_test_file, path_enriched_test_file])
#     kat_run = KatRunner(enriched_sample, path_ref_file, path_output, "test")
#     kat_run.kat_filter(exclude=True)
#     assert kat_run.status == True
#     pass

# def test_kat_hist():
#     enriched_sample = Sequence(technology, [path_enriched_test_file])
#     kat_run = KatRunner(enriched_sample, path_ref_file, path_output, "test")
#     kat_run.kat_hist()
#     assert kat_run.status == True
#     pass

# def test_run_fastp(): 
#     enriched_sample = Sequence(technology, ["/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_1.fastq", "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_2.fastq"])
#     fastp_run = FastPRunner(enriched_sample, path_output, "test_illumina_fastp", report_only=False, dedup=True)
#     fastp_run.run_fastp()
#     assert fastp_run.status == True
#     pass

# def test_run_minimap2(): 
#     enriched_sample = Sequence(technology, ["/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_1.fastq", "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_2.fastq"])
#     minimap2_run = Minimap2Runner(enriched_sample, path_output, "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Saccharomyces_cerevisiae_draft_genome.fasta", "test_illumina_sam")
#     minimap2_run.run_minimap2()
#     assert minimap2_run.status == True
#     pass

# def test_run_samtools_bam(): 
#     samtools_run = SamBamProcessor("/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_illumina_sam.sam", path_output, "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Saccharomyces_cerevisiae_draft_genome.fasta", "test_illumina_bam")
#     samtools_run.run_samtools_bam()
#     assert samtools_run.status == True
#     pass

# def test_run_samtools_fastq(): 
#     samtools_run = SamBamProcessor(bam_test_file, path_output, path_ref_file, "test_output")
#     samtools_run.run_samtools_fastq()
#     assert samtools_run.status == True
#     pass

# def test_run_samtools_import(): 
#     samtools_run = SamBamProcessor(path_enriched_test_file, path_output, path_ref_file, "test_output")
#     samtools_run.run_samtools_import()
#     assert samtools_run.status == True
#     pass

# def test_run_bedtools(): 
#     bedtools_run = SamBamProcessor(bam_test_file, path_output, path_ref_file, "test_output")
#     bedtools_run.run_bedtools(nonzero=True)
#     assert bedtools_run.status == True
#     pass

# def test_run_bam_num_reads(): 
#     bam_run = BamProcessor("/home/ameknas/sequenoscope-1/test/sample_mapped_bam.bam")
#     reads = []
#     for i in bam_run.ref_stats:
#         num_reads = bam_run.ref_stats[i]['num_reads']
#         reads.append(num_reads)
#     print(sum(reads), file=open('test_reads_total.txt', 'a'))
#     assert bam_run.status == True
#     pass

# def test_run_bam(): 
#     bam_run = BamProcessor("/home/ameknas/sequenoscope-1/test/sample_mapped_bam.bam")
#     print(bam_run.ref_stats, file=open('test_output_nano.txt', 'a'))
#     print(bam_run.ref_coverage, file=open('test_output_cov.txt', 'a'))
#     assert bam_run.status == True
#     pass

# def test_run_seq_manifest_with_sum(): 
#     seq_mani_run_nano = SeqManifest("barcode1",
#                                "/home/ameknas/sequenoscope-1/test_SE/sample_mapped_bam.bam", 
#                                "test_out_mani_1",
#                                out_dir=path_output,
#                                fastp_fastq=["/home/ameknas/sequenoscope-1/test_SE/sample_fastp_output.fastp.fastq"],
#                                in_seq_summary= "/home/ameknas/sequenoscope-1/sequenoscope/test_sequences/Nanopore_enriched_control_seq_summary.txt"
#                                )
#     print(seq_mani_run_nano.filtered_reads, file=open('test_output_fastp.txt', 'a'))
#     assert seq_mani_run_nano.status == True
#     pass

# def test_fastq_extractor_sr():
#     enriched_sample = Sequence(technology, ["/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/ERR2984773_1.fastq", "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/ERR2984773_2.fastq"])
#     extractor_run = FastqExtractor(enriched_sample, out_prefix="test_222_reads_sr", out_dir=path_output)
#     extractor_run.extract_paired_reads()
#     assert extractor_run.status == True
#     pass

# def test_fastq_extractor_lr():
#     enriched_sample = Sequence(technology, ["/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/barcode01_fastq_pass_enriched.fastq"])
#     extractor_run = FastqExtractor(enriched_sample, out_prefix="test_reads_lr", out_dir=path_output)
#     extractor_run.extract_single_reads()
#     assert extractor_run.status == True
#     pass


# def test_run_seq_manifest_no_sum(): 
#     seq_mani_run_Ill = SeqManifest("barcode1",
#                                "/home/ameknas/sequenoscope-1/test_SE_no_seq_summary/sample_mapped_bam.bam", 
#                                "test_out_mani_2_Fat_Tits",
#                                out_dir=path_output,
#                                fastp_fastq=["/home/ameknas/sequenoscope-1/test_SE_no_seq_summary/sample_fastp_output.fastp.fastq"],
#                                read_list="/home/ameknas/sequenoscope-1/test_SE_no_seq_summary/sample_read_list.txt",
#                                in_fastq=["/home/ameknas/sequenoscope-1/sequenoscope/test_sequences/barcode01_fastq_pass_enriched.fastq"],
#                                start_time=0,
#                                end_time=100
#                                )
#     #print(seq_mani_run_Ill.filtered_reads.keys(), file=open('test_output_fastp.txt', 'a'))
#     assert seq_mani_run_Ill.status == True
#     pass

# def test_fastq_renamer():
#     enriched_sample = Sequence(technology, ["/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/ERR2984773_1.fastq", "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/ERR2984773_2.fastq"])
#     renamer = FastqPairedEndRenamer(enriched_sample, "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_222_reads_sr.txt", path_output, "test")
#     renamer.rename()
#     renamer.status == True
#     pass

# def test_seq_manifest_summary():
#     seq_mani_run_nano = SeqManifest("barcode1",
#                                "/home/ameknas/sequenoscope-1/test_mani_lr_sum/sample_mapped_bam.bam", 
#                                "test_out_mani_1",
#                                out_dir=path_output,
#                                fastp_fastq=["/home/ameknas/sequenoscope-1/test_mani_lr_sum/sample_fastp_output.fastp.fastq"],
#                                in_seq_summary= "/home/ameknas/sequenoscope-1/sequenoscope/analyze/test_sequences/Nanopore_enriched_control_seq_summary.txt"
#                                )
    
#     kmer_file = GeneralSeqParser("/home/ameknas/sequenoscope-1/test_mani_lr_sum/sample_kmer_analysis_histogram_file.dist_analysis.json", "json")
#     fastp_file = GeneralSeqParser("/home/ameknas/sequenoscope-1/test_mani_lr_sum/sample_fastp_output.json", "json")

#     seq_summary_run = SeqManifestSummary("barcode1",
#                                seq_mani_run_nano.bam_obj, 
#                                "test_out_mani_summary",
#                                out_dir=path_output,
#                                kmer_json_file=kmer_file.parsed_file,
#                                fastp_json_file=fastp_file.parsed_file,
#                                paired=False
#                                )
    
#     seq_summary_run.generate_summary()
#     assert seq_summary_run.status == True
#     pass
    

def test_mash_sketcher(): 
    mash_sketcher = MashSketcher("/home/ameknas/sequenoscope/sequenoscope", "test")
    results_single = mash_sketcher.run_mash_sketch(["/home/ameknas/sequenoscope-1/Sequenoscope/test_sequences/barcode01_fastq_pass_enriched.fastq"])
    print(f"Results for paired-end: {results_single}")

#     # Example usage for paired-end sequencing
#     results_paired = mash_sketcher.run_mash_sketch(["/home/ameknas/sequenoscope-1/Sequenoscope/test_sequences/ERR2984773_1.fastq", "/home/ameknas/sequenoscope-1/Sequenoscope/test_sequences/ERR2984773_2.fastq"])
#     print(f"Results for paired-end: {results_paired}")