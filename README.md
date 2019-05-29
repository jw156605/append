## Installation

- Navigate to desired installation directory
- Download the latest version of AppEnD
- Enter the AppEnD directory
- Download the latest version of the bamtools API
- Build from source

The sequence of commands to do this is:
```
git clone https://github.com/jw156605/append
cd append/
git clone --depth 1 --branch v2.5.1 git://github.com/pezmaster31/bamtools.git
make
```

Note: The above installation process requires that git and cmake be installed on your machine. These can
be obtained from http://git-scm.com/downloads and http://www.cmake.org/download/

## Usage
The AppEnD tool identifies untemplated 3' additions in RNA-seq reads. AppEnD can accurately process many 
different types of sequencing data, provided that the expected position of the untemplated addition within
the read is accurately specified in a parameter file. Sample parameter files are provided for processing 
EnD-Seq data (in which the untemplated additions are present at the beginning of read 1) and A-Seq data 
(single end reads with additions at the end of the read). AppEnD takes sequencing data in the form of an 
indexed BAM file (specified in the parameter file) as input. A crucial part of making AppEnD work is 
using an aligner that performs soft clipping. For example, MapSplice, TopHat, STAR, and bowtie2 (in 
--local mode) can all perform soft clipping. 

To run AppEnD, simply type:
```
./AppEnD parameter_file bam_file
```
