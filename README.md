<p align="center">
  <img src="http://benoukraf-lab.com/wp-content/uploads/2019/05/Nanovarlogo.png" width="200" alt="accessibility text" align='left'>
</p>

# NanoVar - Structural variant caller using low-depth Nanopore sequencing 

Latest version 1.0.1

## Introduction
NanoVar is a neural-network-based genomic structural variant (SV) caller that utilizes low-depth long-read sequencing data generated from Oxford Nanopore Technologies (ONT). It characterizes SVs with high accuracy and speed using only 4x depth genomic sequencing datasets, thereby saving time, sequencing cost, and computational storage, which makes it compatible for large-scale cohort SV-association studies or routine clinical SV investigations.  

### Basic capabilities
* Accurate SV characterization using long sequencing reads (High SV recall and precision in simulation datasets, overall F1 score >0.9)  
* Characterizes many classes of SVs including novel-sequence insertion, genomic-sequence insertion, translocation, deletion, inversion, and tandem duplication  
* Minimal sequencing data required (4x depth whole genome long-read sequencing data)  
* Rapid computational speed (Takes about 3 hours to map and analyze 12 gigabases datasets (4x) using 24 CPU threads)  
* Estimates SV zygosity (Estimates whether an SV is heterozygous or homozygous)  
* Calculates NGS short-read coverage at each SV breakend called by long reads (Work in progress)  


### Software description
* Written in Shell and Python  
* Utilizes HS-BLASTN (Chen et al., 2015) for long-read mapping against a reference genome  
* Employs Keras with TensorFlow backend for neural network inferencing for SV confidence calling  
* Creates a Python virtual environment for the execution of all Python processes including Keras  


## Getting Started

### 1. Requirements and Prerequisites

#### Operating system requirements: 
* Linux (x86_64 architecture, tested in Ubuntu 16.04 and Ubuntu 14.04)  
* Mac (Work in progress)  

#### Prerequisites:
1. Python 2 (version 2.7.12 or higher)  
2. Samtools (version 1.4 or higher)  
3. Bedtools (version 2.26 or higher)  

Optional: Bowtie2 (version 2.3.4 or higher, required only if short-reads will be incorporated, otherwise, do "./configure --disable-bowtie2" during installation)  


### 2. Installation
#### Clone git repository:
```
git clone https://github.com/benoukraflab/NanoVar.git 
cd NanoVar 
./configure
make && make check
sudo make install # or add the executable "nanovar" to PATH
```
Or 

#### Download source code:
Download tarball or zipped source code from [Releases](https://github.com/benoukraflab/NanoVar/releases)
```
tar zxvf NanoVar-x.x.tar.gz # or unzip NanoVar-x.x.zip
cd NanoVar-x.x
./configure
make && make check
sudo make install # or add the executable "nanovar" to PATH
```

### 3. Quick run

#### Run with long reads only
```
nanovar [Options] -t 24 -r hg38.fa -l longread.fa -o ./output/ -f hg38 
```
#### Run with long reads and NGS paired-end short reads
```
nanovar [Options] -t 24 -r hg38.fa -l longread.fa -s1 shortread_mate1.fq -s2 shortread_mate2.fq -o ./output/ -f hg38
```

| Parameter | Argument | Comment |
| :--- | :--- | :--- |
| `-t` | num_threads | Indicate number of CPU threads to use |
| `-r` | reference.fa | Input reference genome in FASTA format |
| `-l` | longread.fa | Input long-read FASTA/FASTQ file |
| `-o` | output_directory | Define output directory |
| `-f` | gap_file | Choose built-in gap BED file to exclude gap regions in the reference genome. Built-in gap files include: hg19, hg38, mm9, and mm10 |
| `-s1` | shortread_mate1.fq | NGS short-read paried-end mate 1 FASTA/FASTQ file |
| `-s2` | shortread_mate2.fq | NGS short-read paried-end mate 2 FASTA/FASTQ file |

## Documentation
See [Wiki](https://github.com/benoukraflab/NanoVar/wiki) for more information.

## Versioning
See [Releases](https://github.com/benoukraflab/NanoVar/releases)

For development releases, please visit [nanovar-dev](https://github.com/cytham/nanovar-dev)

## Citation

## Authors

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)
* **Roberto Tirado Magallanes** - [rtmag](https://github.com/rtmag)
* **Touati Benoukraf** - [benoukraflab](https://github.com/benoukraflab)

## License

This project is licensed under GNU General Public License - see [COPYING](https://github.com/benoukraflab/NanoVar/blob/master/COPYING) for details.

## Simulation datasets
SV-simulated datasets used for evaluating SV calling accuracy can be downloaded [here](https://doi.org/10.5281/zenodo.2599376).
