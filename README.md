[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/sequenoscope/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/sequenoscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sequenoscope/badges/downloads.svg)](https://anaconda.org/bioconda/sequenoscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sequenoscope/badges/license.svg)](https://anaconda.org/bioconda/sequenoscope)


## Sequenoscope

<img height="150" width="400" alt="logo11" src="https://user-images.githubusercontent.com/93303799/225326096-6c9de0f1-9ac0-46a4-914f-a3db51e7e97a.png">

A tool for analyzing sequencing run outputs primarily from adaptive sampling experiments and Oxford Nanopore Technology sequencers.

## Contents

- [Introduction](#introduction)
- [Dependencies](#Dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Use-case Example](#Use-case-example)
- [Quick Start](#quick-start)
  -  [Analyze](#analyze-module)
  -  [Filter_ONT](#filter_ONT-module)
  -  [Plot](#plot-module)
- [Outputs](#Outputs)
  -  [Analyze](#analyze-module-outputs)
      -  [Manifest file report structure](#sample-manifest-report-format)
      -  [Manifest summary file report structure](#sample-manifest-summary-report-format)
  -  [Filter_ONT](#filter_ONT-module-outputs)
  -  [Plot](#plot-module-outputs)
- [Benchmarks](#Benchmarks)
- [FAQ](#faq)
- [Citation](#citation)
- [Legal](#legal)
- [Contact](#contact)

## Introduction

Analyzing and interpreting sequencing data is a fundamental task in bioinformatics, and with the advent of ONT adaptive-sampling sequencing, specialized tools are needed to visualize and assess sequencing runs. The unique data output and complex dynamics of adaptive sampling present challenges in effectively visualizing and assessing these sequencing runs in terms of key parameters, necessitating tailored analytical approaches and visual analytics. To address these challenges, we have developed a comprehensive bioinformatics pipeline consisting of three modules:  [analyze](#analyze-module), [plot](#plot-module), and [filter_ONT](#filter_ONT-module). Our accessible pipeline aims to provide researchers with a fast and intuitive workflow for easily processing and analyzing sequencing data especially from adaptive sequencing runs, enabling them to gain interpretable insights into their datasets with minimal upfront efforts.

The [analyze](#analyze-module) module serves as the core component of our pipeline. It takes an input FASTQ file, a reference FASTA file, and an optional sequencing summary file from ONT sequencers. Leveraging tools such as `fastp`, `minimap2`, `pysam`, and `mash`, this module performs a series of essential tasks. It filters the input FASTQ file, maps it to the reference FASTA file, and generates a **sequence manifest file** and **summary sequence manifest file**. These files include key sequencing [statistics](#sample-manifest-report-format) such as read length, read quality (Q score), mapping efficiency, and coverage depth. For an in-depth explanation of all statistics provided, please refer to the [report format section](#file-report-section) below.

The [plot](#plot-module) module complements the analysis performed by the "analyze" module. It takes as input both a "test" and "control" directory containing __manifest__ and __manifest summary files__ generated by the "analyze" module. With these files, the [plot](#plot-module) module generates visualizations that aid in the interpretation and visualization of the sequencing data. Due to its focus on comparative analysis, a quality control (QC) run that does not include distinct "test" or "control" samples would not be applicable for analysis with this module.

The [filter_ONT](#filter_ONT-module) module is designed for for raw reads filtering and subsetting. This module leverages a sequencing summary file to allow researchers to precisely filter reads based on customized criteria, including channel, duration, start time, [Q score](https://community.nanoporetech.com/technical_documents/data-analysis/v/datd_5000_v1_revr_22aug2016/basecalling-considerations), and length.

Our bioinformatics pipeline offers a powerful tool for researchers working with ONT sequencing data. Whether you are exploring metagenomics sample composition, investigating adaptive sampling for your project, or conducting a comparative analysis of different methods in your lab, our pipeline can streamline your analyses and provide valuable insights into your genomic datasets using visual aids and easy to understand outputs.

## Dependencies
- Python: `>=3.7.12, <4`
- fastp: `>=0.22.0`
- mash: `>=2.3`
- minimap2: `>=2.26`
- seqtk: `>=1.4`
- samtools: `>=1.6`

## Python Packages
- pysam: `>=0.16.0`
- plotly: `>=5.16.1`

## Installation

## Option 1: As a conda package (Recomended)

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n sequenoscope

## Option 2: As a PyPI package

**Coming soon**

Install using pip:

        pip install sequenoscope

## Option 3: Install from source

**Coming soon**

If you wish to install sequenoscope from source, please first ensure these dependencies are installed and configured on your system:
`python>=3.7.12,<4`
`fastp >=0.22.0`
`mash >=2.3`
`minimap2 >=2.26`
`seqtk >=1.4`
`samtools >=1.6`
`pysam >=0.16.0`
`plotly >=5.16.1`

Install the latest commit from the master branch directly from Github:

        pip install git+https://github.com/phac-nml/sequenoscope.git

## Usage
If you run ``sequenoscope``, you should see the following usage statement:

        Usage: sequenoscope <command> <required arguments>
        
        To get full help for a command use one of:
        sequenoscope <command> -h
        sequenoscope <command> --help
        
        
        Available commands:
        
        analyze     map reads to a target and produce a report with sequencing statistics
        plot        generate plots based on directories with seq manifest files
        filter_ONT  filter reads from a FASTQ file based on a sequencing summary file


## Use-case example
To demonstrate the practical application of our pipeline, consider a scenario where a researcher conducts adaptive sampling using an ONT sequencer. In this example, the researcher divides the sequencer channels into two sets: one half for adaptive sampling enrichment and the other half for regular sequencing as a control.

Utilizing our [filter_ONT](#filter_ONT-module) module, the researcher can create two distinct sets of FASTQ files, each representing the minimum and maximum channels of the sequencing data.

These files are then processed separately through our [analyze](#analyze-module) module, generating two datasets – one for the test (adaptive sampling) and one for the control (regular sequencing).

Finally, by employing the [plot](#plot-module) module, the researcher can visually assess the effectiveness of the adaptive sampling enrichment in their experiment. This practical illustration highlights how our pipeline facilitates comprehensive data processing and analysis, enhancing the researcher's ability to draw meaningful conclusions from their ONT sequencing data.

## Quick start

## analyze module

To quickly get started with the `analyze` module, follow the instructions below:

1. Ensure that you have the necessary input files and reference database prepared:
   - **Input FASTQ files:** Provide the path to the FASTQ files you want to process using the `--input_FASTQ` option.
   - **Reference database:** Specify the path to the reference database in FASTA format using the `--input_reference` option.

2. Choose an output directory for the results:
   - Specify the output directory using the `--output` option. This is a required parameter.

3. Run the command with the basic required options:

        sequenoscope analyze --input_fastq <file.fq> --input_reference <ref.FASTA> -o <output_directory> -seq_type <sr>

This command will initiate the analysis process using default settings. The FASTQ files will be processed, and the results will be saved in the specified output directory.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope analyze -h` or `sequenoscope analyze --help`.

Remember to replace `<file.fq>` with the actual path to your FASTQ file, `<ref.FASTA>` with the path to your reference database, `<output_directory>` with the desired location for the output files and `<sr>` with your sequencing type (SE for single-end and PE for paired-end).

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

## filter_ONT module

To quickly get started with the `filter_ONT` module, follow the instructions below:

1. Ensure that you have the necessary input files prepared:
   - **Input FASTQ files:** Provide the path to the adaptive sequencing FASTQ files you want to process using the `--input_FASTQ` option.
   - **ONT sequencing summary file:** Specify the path to the ONT sequencing summary file using the `--input_summary` option.

2. Choose an output file and directory for the filtered reads:
   - Specify the output file path and directory using the `--output` option. This is a required parameter.

3. Set the desired filtering criteria:
   - You can apply various filters to the reads based on the following options:
     - Classification: Use the `-cls` or `--classification` option to designate the adaptive-sampling sequencing decision classification. Valid options are `'unblocked'`, `'stop_receiving'`, or `'no_decision'`.
     - Channel/Pore number: Set the minimum and maximum channel/pore number for filtering using the `-min_ch` and `-max_ch` options.
     - Duration: Define the minimum and maximum duration of the sequencing run in seconds using the `-min_dur` and `-max_dur` options.
     - Start time: Specify the minimum and maximum start time of the sequencing run in seconds using the `-min_start` and `-max_start` options.
     - Q score: Determine the minimum and maximum q score for filtering using the `-min_q` and `-max_q` options.
     - Read length: Set the minimum and maximum read length for filtering using the `-min_len` and `-max_len` options.

4. Run the command with the basic required options:

        sequenoscope filter_ONT --input_fastq <file.fq> --input_summary <seq_summary.txt> -o <output.FASTQ>

This command will initiate the filtering process based on the specified criteria and save the filtered reads to the output file.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope filter_ONT -h` or `sequenoscope filter_ONT --help`.

Remember to replace `<file.fq>` with the actual path to your ONT sequencing FASTQ file, `<seq_summary.txt>` with the path to your ONT sequencing summary file, and `<output.FASTQ>` with the desired path and filename for the filtered reads.

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

## plot module
To quickly get started with the `plot` module, follow the instructions below:

1. **Required Paths:** Ensure you have designated the necessary directories:
- **Test Directory:** Provide the path to the test directory that contains the __seq manifest__ files from the analyze module. `-T` or `--test_dir` `<test_dir_path>`
- **Control Directory:** Specify the path to the control directory that contains the __seq manifest__ files from the analyze module. `-C` or `--control_dir` `<control_dir_path>`
- **Output Directory:** Choose an output directory for the plots. `-o` or `--output_dir` `<out_path>`

2. **Plotting Options:** Customize your plots with various options:
- **Output Prefix:** You can add a prefix before plot names with the `--output_prefix` option. `-op` or `--output_prefix` `<OUTPUT_PREFIX>`. *Default is 'sample'*.
- **Comparison Metric:** Select a parameter for the box plot and single ratio bar chart using the `--comparison_metric` option. *Default parameter is taxon_%_covered_bases*.
- **Single Charts:** Generate charts for data based on selected comparison metric using the `--single_charts` option. `{TRUE, FALSE}`. *Default designation is False*
- **Adaptive Sampling:** Generate decision bar charts for adaptive sampling if utilized during sequencing using the `-AS` option. `{TRUE, FALSE}`. *Default designation is False*
- **Violin Data Fraction:** Set a fraction of the data to use for the violin plot. `-VP` or `--violin_data_percent` `<0.1 - 1>`. *Default fraction is 0.1*
- **Time Bin Unit:** Designate a time bin used for decision bar charts. `-bin` or `--time_bin_unit` `{seconds,minutes,5m,15m,hours}`. *Default bin is minutes*
- **Taxon Legend:** Generate a legend for the source file taxon covered bar chart. `-legend` or `--taxon_chart_legend` `{TRUE, FALSE}`. *Default designation is False*

3. **Run the Command:** With the basic required options:

        sequenoscope plot --test_dir <test_dir_path> --control_dir <control_dir_path> --output_dir <out_path>

Use the `--force flag` if you wish to force an overwrite of an existing results directory.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope plot -h` or `sequenoscope plot --help`.

Remember to replace `<test_dir_path>`, `<control_dir_path>`, and `<out_path>` with the actual paths for your directories.

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

## Outputs

## analyze module outputs

| File | Description |
|------|-------------|
| `<prefix>_fastp_output.fastq` | The output FASTQ file after processing with `fastp`. It includes filtered and trimmed sequencing reads. |
| `<prefix>_fastp_output.html` | An HTML report generated by `fastp` summarizing the filtering and quality control results. |
| `<prefix>_fastp_output.json` | A JSON formatted report with detailed `fastp` quality control statistics. |
| `<prefix>_manifest.txt` | A sequence manifest file containing various sequencing statistics post-analysis. |
| `<prefix>_manifest_summary.txt` | A summary of the sequence manifest with key statistics for a quick overview. |
| `<prefix>_mapped.bam` | The BAM file output from `minimap2`, containing aligned sequences to the reference FASTA. |
| `<prefix>_mapped.bam.bai` | An index file for the BAM file to enable quick read access. |
| `<prefix>_mapped_fastq.fastq` | The FASTQ file containing reads that have been mapped to the reference. |
| `<prefix>_mapped.sam` | The SAM file equivalent of the BAM file, containing human-readable alignment data. |
| `<prefix>_mash.hash.msh` | A MASH sketch file used for rapid genome distance estimation. |
| `<prefix>_read_list.txt` | A text file list of reads, potentially used for further downstream analysis. |

Note: Replace `<prefix>` with the user-specified prefix that precedes all output filenames.

### sample manifest report format

| Column ID | Description |
|-----------|-------------|
| `sample_id` | Identifier for the sample to which the read belongs. |
| `read_id` | Unique identifier for the sequencing read. |
| `read_len` | Length of the sequencing read in base pairs. |
| `read_qscore` | Quality score of the sequencing read. |
| `channel` | The channel on the sequencing device from which the read was recorded. |
| `start_time` | Time when the sequencing of the read started. |
| `end_time` | Time when the sequencing of the read ended. |
| `decision` | Indicates the final decision on the sequencing read. Decisions are categorized into three main types: `stop_receiving` (the sequencing is allowed to continue, represented by `signal_positive`), `unblocked` (the read is ejected from sequencing, indicated by `data_service_unblock_mux_change`), and `no_decision` (no definitive action was taken, denoted by either `signal_negative` or `unblock_mux_change`). Each term explains the action taken or not taken based on the read's signal detection and processing status. |
| `fastp_status` | Indicates whether the read passed the filtering and trimming process by `fastp`. |
| `is_mapped` | Indicates whether the read is mapped to any sequence in the provided multi-sequence FASTA reference file (`TRUE` if mapped, also see note 1 below). |
| `is_uniq` | Indicates whether the read is unique within the sample manifest file (`TRUE` if unique, also see note 2 below). |
| `contig_id` | Identifier for the contig to which the read is mapped, if applicable. |

**Notes:**
1. `is_mapped` refers to whether or not a read is mapped to any sequence in the multi-sequence FASTA reference file provided by the user. If true, the `contig_id` is provided.
2. `is_uniq` refers to whether or not a read is unique throughout the sample manifest file. In ONT sequencing, a read may be processed multiple times if the decision is labelled as `signal_negative` or `No_decision` before a final decision is made on whether to allow the read to continue sequencing or not.

### sample manifest summary report format

| Column ID | Description |
|-----------|-------------|
| `sample_id` | Identifier for the sample. |
| `est_genome_size` | Estimated size of the genome. |
| `est_coverage` | Estimated coverage of the genome. |
| `total_bases` | Total number of bases in the sample. |
| `total_fastp_bases` | Total number of bases after processing with `fastp`. |
| `mean_read_length` | Mean read length of the sequencing reads. |
| `taxon_id` | Identifier for the taxon. Obtained from the user-provided FASTA file. |
| `taxon_length` | Length of the taxon's genome. |
| `taxon_mean_coverage` | Mean coverage across the taxon's genome. |
| `taxon_covered_bases_<prefix>X` | Number of bases in the taxon's genome covered at user-specified coverage threshold. |
| `taxon_%_covered_bases` | Percentage of the taxon's genome that is covered by reads at the user-specified coverage threshold . |
| `total_taxon_mapped_bases` | Total number of bases mapped to the taxon. |
| `taxon_mean_read_length` | Mean read length of the reads mapped to the taxon. |


Note: Replace `<prefix>` with the user-specified threshold coverage.

## filter_ONT module outputs

| File | Description |
|------|-------------|
| `<user_prefix>_filtered_fastq_subset.fastq` | The subset of FASTQ reads that have been filtered based on the user-defined criteria within the `filter_ONT` module. |
| `<user_prefix>_read_id_list.csv` | A CSV file containing the list of read identifiers that correspond to the filtered subset. This may be used for further reference or analysis. |

Note: Replace `<prefix>` with the user-specified prefix that precedes all output filenames.

## plot module outputs

| File | Description | Triggered by Command |
|------|-------------|----------------------|
| `<prefix>_ratio_bar_chart.html` | An HTML file containing a bar chart that displays the ratio statistics of the manifest summary file. | Default behavior |
| `<prefix>_source_file_taxon_covered_bar_chart.html` | An HTML file containing a bar chart displaying the coverage of taxa in the source files. | Default behavior and `--taxon_chart_legend` specifying the inclusion of a legend|
| `<prefix>_stat_results.csv` | A CSV file with statistical results of the analysis, such as taxa coverage percentages. | Default behavior |
| `<prefix>_cumulative_decision_bar_chart.html` | An HTML file containing a bar chart with cumulative decision metrics over time for either test or control datasets. | adaptive sampling enabled (`-AS`) and time-bin specified (`--time_bin_unit`)|
| `<prefix>_independent_decision_bar_chart.html` | An HTML file containing a bar chart with independent decision metrics over time for either test or control datasets. | adaptive sampling enabled (`-AS`) and time-bin specified (`--time_bin_unit`) |
| `read_len_<prefix>_violin_comparison_plot.html` | An HTML file containing a violin plot comparing log-transformed data between the test and control datasets. | Default behavior and `--violin_data_percent` specifying the fraction of data to plot |
| `read_qscore_<prefix>_violin_comparison_plot.html` | An HTML file containing a violin plot comparing q-score distributions between test and control datasets. | Default behavior and `--violin_data_percent` specifying the fraction of data to plot |
| `<prefix>_box_plot.html` | Generate a box plot comparing a specific parameter from test and control files. |  `--comparison_metric` specified with `--single_charts` enabled |
| `<prefix>_single_ratio_bar_chart.html` | Generate a single bar chart comparing a specific parameter from test and control files. |  `--comparison_metric` specified with `--single_charts` enabled |

Note: Replace `<prefix>` with the user-specified prefix that precedes all output filenames from the `plot` module. This prefix is set with the `--output_prefix` option when running the command.

Note: For the adaptive sampling plots specified with `-AS` command, there will be 2 files, test and control, for each type of bar chart, independent and cumulative.

## Benchmarks

**Coming soon**

## FAQ

**Coming soon**

## Citation

**Coming soon**

## Legal

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.


## Contact

**Abdallah Meknas**: abdallah.meknas@phac-aspc.gc.ca

