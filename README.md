[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/sequenoscope/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/sequenoscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sequenoscope/badges/downloads.svg)](https://anaconda.org/bioconda/sequenoscope)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sequenoscope/badges/license.svg)](https://anaconda.org/bioconda/sequenoscope)


## Sequenoscope

<img height="150" width="400" alt="logo11" src="https://user-images.githubusercontent.com/93303799/225326096-6c9de0f1-9ac0-46a4-914f-a3db51e7e97a.png">

A tool for analyzing sequencing output. 

## Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [FAQ](#faq)
- [Citation](#citation)
- [Legal](#legal)
- [Contact](#contact)

## Introduction

Analyzing and interpreting sequencing data is a fundamental task in bioinformatics, especially when studying adaptive sampling MinION sequencing data. To address this challenge, we have developed a comprehensive bioinformatics pipeline consisting of three modules: "analyze," "plot," and "filter_ONT." Our pipeline aims to provide researchers with a streamlined and efficient workflow for processing and analyzing sequencing data, enabling them to gain valuable insights into their datasets.

The "analyze" module serves as the core component of our pipeline. It takes an input FASTQ file, a reference FASTA file, and an optional sequencing summary file from ONT sequencers. Leveraging tools such as fastp, minimap2, pysam, and mash, this module performs a series of essential tasks. It filters the input FASTQ file, maps it to the reference FASTA file, and generates sequence manifest files and summary sequence manifest files, which contain various sequencing statistics.

The "plot" module complements the analysis performed by the "analyze" module. It takes as input both a test and control directory containing manifest and manifest summary files generated by the "analyze" module. With these files, the "plot" module generates visualizations that aid in the interpretation and visualization of the sequencing data.

To further enhance the versatility of our pipeline, we have developed the "filter_ONT" module. This module enables users to filter reads from a FASTQ file based on a sequencing summary file. Researchers can customize their filtering criteria by specifying parameters such as channel, duration, start time, q score, and length.

To demonstrate the practical application of our pipeline, consider a scenario where a researcher conducts adaptive sampling using a MinION sequencer. In this example, the researcher divides the sequencer pores into two sets: one half for adaptive sampling enrichment and the other half for regular sequencing as a control. Utilizing our "filter_ONT" module, the researcher can create two distinct sets of FASTQ files, each representing the minimum and maximum channels of the sequencing data. These files are then processed separately through our "analyze" module, generating two datasets – one for the test (adaptive sampling) and one for the control (regular sequencing). Finally, by employing the "plot" module, the researcher can visually assess the effectiveness of the adaptive sampling enrichment in their experiment. This practical illustration highlights how our pipeline facilitates comprehensive data processing and analysis, enhancing the researcher's ability to draw meaningful conclusions from their MinION sequencing data.

By providing a comprehensive suite of modules, our bioinformatics pipeline offers a powerful tool for researchers working with adaptive sampling MinION sequencing data. Whether you are exploring metagenomics sample composition, investigating adaptive sampling for your project, or conducting a comparative analysis of different methods in your lab, our pipeline can streamline your analyses and provide valuable insights into your genomic datasets.

## Installation

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n sequenoscope

**Coming soon**

Install using pip:

        pip install sequenoscope

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

Install the latest master branch version directly from Github:

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

Quick start
=====
To quickly get started with the `analyze` module, follow the instructions below:

1. Ensure that you have the necessary input files and reference database prepared:
   - **Input FASTQ files:** Provide the path to the FASTQ files you want to process using the `--input_FASTQ` option.
   - **Reference database:** Specify the path to the reference database using the `--input_reference` option.

2. Choose an output directory for the results:
   - Specify the output directory using the `--output` option. This is a required parameter.

3. Run the command with the basic required options:

        sequenoscope analyze --input_FASTQ <file.fq> --input_reference <ref.FASTA> -o <output_directory> -seq_type <sr>

This command will initiate the analysis process using default settings. The FASTQ files will be processed, and the results will be saved in the specified output directory.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope analyze -h` or `sequenoscope analyze --help`.

Remember to replace `<file.fq>` with the actual path to your FASTQ file, `<ref.FASTA>` with the path to your reference database, `<output_directory>` with the desired location for the output files and `<sr>` with your sequencing type (SE for single-end and PE for paired-end).

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

-------

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

        sequenoscope filter_ONT --input_FASTQ <file.fq> --input_summary <seq_summary.txt> -o <output.FASTQ>

This command will initiate the filtering process based on the specified criteria and save the filtered reads to the output file.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope filter_ONT -h` or `sequenoscope filter_ONT --help`.

Remember to replace `<file.fq>` with the actual path to your adaptive sequencing FASTQ file, `<seq_summary.txt>` with the path to your ONT sequencing summary file, and `<output.FASTQ>` with the desired path and filename for the filtered reads.

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

-------------
To quickly get started with the `plot` module, follow the instructions below:

1. **Required Paths:** Ensure you have designated the necessary directories:
- **Test Directory:** Provide the path to the test directory that contains the seq manifest files from the analyze module. `-T` or `--test_dir` `<test_dir_path>`
- **Control Directory:** Specify the path to the control directory that contains the seq manifest files from the analyze module. `-C` or `--control_dir` `<control_dir_path>`
- **Output Directory:** Choose an output directory for the plots. `-o` or `--output_dir` `<out_path>`

2. **Plotting Options:** Customize your plots with various options:
- **Output Prefix:** You can add a prefix before plot names with the `--output_prefix` option. `-op` or `--output_prefix` `<OUTPUT_PREFIX>`. *Default is 'sample'*.
- **Comparison Metric:** Select a parameter for the box plot and single ratio bar chart using the `--comparison_metric` option. *Default parameter is taxon_%_covered_bases*.
- **Single Charts:** Generate charts for data based on selected comparison metric using the `--single_charts` option. *Default designation is False*
- **Adaptive Sampling:** Generate decision bar charts for adaptive sampling if utilized during sequencing using the `-AS` option. *Default designation is False*
- **Violin Data Fraction:** Set a fraction of the data to use for the violin plot. `-VP` or `--violin_data_percent` `<0.1 - 1>`. *Default fraction is 0.1*
- **Time Bin Unit:** Designate a time bin used for decision bar charts. `-bin` or `--time_bin_unit` `{seconds,minutes,5m,15m,hours}`. *Default bin is minutes*

3. **Run the Command:** With the basic required options:

        sequenoscope plot --test_dir <test_dir_path> --control_dir <control_dir_path> --output_dir <out_path>

Use the `--force flag` if you wish to force an overwrite of an existing results directory.

Please note that this is a simplified quick start guide, and additional options are available for advanced usage. For more detailed information on available options, you can run `sequenoscope plot -h` or `sequenoscope plot --help`.

Remember to replace `<test_dir_path>`, `<control_dir_path>`, and `<out_path>` with the actual paths for your directories.

If you encounter any issues or need further assistance, refer to the full documentation or consult the available resources for troubleshooting.

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

