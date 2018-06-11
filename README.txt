program     methylCtools
purpose     analysis of whole-genome bisulfite sequencing data

version     1.0.0
date        10 june 2018
author      volker hovestadt
institute   developed at the german cancer research center (dkfz), 2011-2015
contact     methylctools@hovestadt.bio


methylCtools ..
 .. is a set of tools for the analysis of whole-genome bisulfite sequencing data for studying DNA methylation [1].
 .. uses a similar alignment approach as bismark [2] (and others before), but improves on the handling of large amounts of data by providing enhanced speed, scalability and avoiding unnecessary I/O.
 
 .. builds upon the popular short-read aligner BWA [3], integrates well with the most of the common/popular tools of next-generation sequencing analysis. 
 .. implements a novel approach for the detection of SNVs within CpG/non-CpG sites that prevents erroneous methylation calling. check the bcall.py file for more details.

 .. does not provide a one-touch workflow, but requires the user to perform and understand the individual steps.
 .. has been designed for use on large server infrastructures.
 
 

INSTALL:

make sure pysam is working within python and make methylCtools executable.
methylCtools has been tested on linux and mac systems.

dependencies:
python 2    version tested: 2.7.2               http://www.python.org/
pysam       version tested: 0.5                 http://code.google.com/p/pysam/                 (version 0.6 has not been tested yet, but should work in principle)
bwa 0.6     version tested: 0.6.1-r104          http://bio-bwa.sourceforge.net/                 [3] (version >0.6 required for human genome because of reference size limit in previous versions)
samtools    version tested: 0.1.18 (r982:295)   http://samtools.sourceforge.net/                [4]
tabix       version tested: 0.2.5 (r964)        http://sourceforge.net/projects/samtools/files/tabix/

suggestions:
FastqC      version tested: 0.10.0              http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
mbuffer     version tested: 20110724            http://www.maier-komor.de/mbuffer.html
Picard      version tested: 1.61(1094)          http://picard.sourceforge.net/



SYNOPSIS:

- performed once (generation of indexes):

1a)  methylCtools fapos reference.fa - | bgzip > reference.pos.gz
1b)  tabix -s 1 -b 2 -e 2 positions.pos.gz

2a)  methylCtools faconv reference.fa reference.conv.fa
2b)  bwa index -a bwtsw reference.conv.fa

- perfomed for every sample (processing of sequencing reads):

3a)  methylCtools fqconv -1 reads1.fq reads1.conv.fq
3b)  methylCtools fqconv -2 reads2.fq reads2.conv.fq

4)   bwa mem -M reference.conv.fa reads1.conv.fq reads2.conv.fq | samtools view -Sb - > reads.conv.bam

5a)  methylCtools bconv reads.conv.bam - | samtools sort - reads
5b)  samtools index reads.bam

6a)  methylCtools bcall reference.pos.gz sample.bam - | bgzip > sample.call.gz
6b)  tabix -s 1 -b 2 -e 2 sample.call.gz



COMMENTS:

- all tools in methylCtools display help messages that explain different options and input/output files (--help).

- an overview of the different processing steps is provided in FLOWCHART.pdf

- most tools support reading from stdin and writing to stdout. simply define "-" instead of the filename.

- if bwa-mem is used (recommended), mark shorter split reads as secondary (-M) for methylCtools bconv compatibility.

- using pipes significantly speeds up processing time and reduces I/O. for example, steps 3 to 5a can be run in a single pipe (define both -1 and -2 in methylCtools fqconv to create an interleaved fastq file, then use -p in bwa mem, then samtools view, methylCtools bconv, and samtools sort). use mbuffer between tools to minimize waits.

- if you are planning to use picard (for example to remove PCR duplicates), make sure to have the entries in the original fasta file sorted alphabetically (e.g. chr1, chr10, chr11, ..., chrM, chrX, chrY).

- please do not hesitate to send comments, bug reports or feature requests.


CITATION:

V Hovestadt, DTW Jones, S Picelli, W Wang, M Kool, PA Northcott, et al. Decoding the regulatory landscape of medulloblastoma using DNA methylation sequencing. Nature 510 (7506), 537-541 135, 2014



REFERENCES:

[1]  Lister, R., & Ecker, J. R. (2009). Finding the fifth base: genome-wide sequencing of cytosine methylation. Genome Research, 19(6), 959–966. doi:10.1101/gr.083451.108
[2]  Krueger, F., & Andrews, S. R. (2011). Bismark: A flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics (Oxford, England). doi:10.1093/bioinformatics/btr167
[3]  Li, H., & Durbin, R. (2010). Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics (Oxford, England), 26(5), 589–595. doi:10.1093/bioinformatics/btp698
[4]  Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. doi:10.1093/bioinformatics/btp352



LICENSE:

MIT License

Copyright (c) 2011-2018 Volker Hovestadt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

