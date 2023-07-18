[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/sequenoscope/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/sequenoscope)
[![Conda](https://img.shields.io/conda/dn/bioconda/profile_dists?color=green)](https://anaconda.org/bioconda/sequenoscope)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/profile_dists)](https://www.apache.org/licenses/LICENSE-2.0)


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

The "analyze" module serves as the core component of our pipeline. It takes an input fastq file, a reference fastq file, and an optional sequencing summary file from ONT sequencers. Leveraging tools such as fastp, minimap2, pysam, and kat, this module performs a series of essential tasks. It filters the input fastq file, maps it to the reference fastq file, and generates sequence manifest files and summary sequence manifest files, which contain various sequencing statistics. Additionally, it produces informative kmer plots using the kat program.

The "plot" module complements the analysis performed by the "analyze" module. It takes as input both test and control sequence manifest or a sequence manifest summary files generated by the "analyze" module. With these files, the "plot" module generates visualizations that aid in the interpretation and visualization of the sequencing data.

To further enhance the versatility of our pipeline, we have developed the "filter_ONT" module. This module enables users to filter reads from a fastq file based on a sequencing summary file. Researchers can customize their filtering criteria by specifying parameters such as channel, duration, start time, q score, and length.

By providing a comprehensive suite of modules, our bioinformatics pipeline offers a powerful tool for researchers working with adaptive sampling MinION sequencing data. Whether you are exploring metagenomics sample composition, investigating adaptive sampling for your project, or conducting a comparative analysis of different methods in your lab, our pipeline can streamline your analyses and provide valuable insights into your genomic datasets.

## Installation

Install the latest released version from conda:

        conda create -c bioconda -c conda-forge -n sequenoscope

**Coming soon**

Install using pip:

        pip install sequenoscope

**Coming soon**

Install the latest master branch version directly from Github:

        pip install git+https://github.com/phac-nml/sequenoscope.git



## Usage
**Coming soon**

Quick start
=====
**Coming soon**
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

