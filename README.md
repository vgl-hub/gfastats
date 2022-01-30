# gfastats
A single fast and exhaustive tool for **summary statistics** and simultaneous \*fast\* (fasta, fastq, gfa [.gz]) genome assembly file **manipulation**.
**gfastats** also allows seamless fasta<>fastq<>gfa[.gz] conversion.

Metrics include:
- scaffold, contig and gap size
- number of scaffolds, contigs and gaps
- total length of scaffolds, contigs and gaps
- scaffold, contig, gap N50 and statistics (full N\*/NG\* statistics with the `--nstar-report` flag)
- area under the curve (AuN/AuNG) values for scaffolds, contigs and gaps
- average scaffold, contig, gap size
- largest scaffold, contig and gap
- base composition and GC content
- soft-masked base counts (lower case bases)

Metrics for each scaffold/contig can be generated with the `--seq-report` flag.

`Bed` coordinates ans sizes of scaffolds, contigs and gaps can be outputted with the options `--out-coord` and `--out-size`. By default, `--out-coord` produces a full representation of the assembly in `agp` format.

Additionally, input can be filtered using scaffold lists or `bed` coordinate files with the options `--include-bed` and `--exclude-bed`.

Importantly, the filtered input can be outputted in any \*fast\* (fasta, fastq, gfa [.gz]) format.

## Installation
Either download one of the releases or `git clone https://github.com/vgl-hub/gfastats.git` and `make`.

## Usage
`gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]`
To check out all option use `gfastats -h`.

## Description
Please refer to **gfastats** paper for a complete description. Briefly, **gfastats** reads and stores any fasta<>fastq<>gfa[.gz] in gfa format. **gfastats** then builds a bidirected graph representation of the assembly using adjaciency lists, where each node is a segment and each edge is a gap (see figure below). The original sequence can be directly manipulated from the graph. Finally, walking the graph allows to generate different kinds of outputs, including manipulated assemblies and feature coordinates.

## How to cite
If you use **gfastats** in your work please cite:

 
