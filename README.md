# gfastats
The swiss army knife for genome assembly.

**gfastats** is a single fast and exhaustive tool for **summary statistics** and simultaneous \*fa\* (fasta, fastq, gfa [.gz]) genome assembly file **manipulation**.
**gfastats** also allows seamless fasta<>fastq<>gfa[.gz] conversion.

Typical fast\* metrics include:
- scaffold, contig and gap size
- number of scaffolds, contigs and gaps
- total length of scaffolds, contigs and gaps
- scaffold, contig, gap N50 and statistics (full N\*/NG\* statistics with the `--nstar-report` flag)
- area under the curve (AuN/AuNG) values for scaffolds, contigs and gaps
- average scaffold, contig, gap size
- largest scaffold, contig and gap
- base composition and GC content
- soft-masked base counts (lower case bases)

Typical gfa metrics include:
- Number of nodes and edges

Metrics for each scaffold/contig can be generated with the `--seq-report` flag.

`Bed` coordinates and sizes of scaffolds, contigs and gaps can be outputted with the options `--out-coord` and `--out-size`. By default, `--out-coord` produces a full representation of the assembly in `agp` format.

Additionally, input can be filtered using scaffold lists or `bed` coordinate files with the options `--include-bed` and `--exclude-bed`.

Importantly, the filtered input can be outputted in any \*fa\* (fasta, fastq, gfa [.gz]) format.

## Installation
Either download one of the releases or `git clone https://github.com/vgl-hub/gfastats.git` and `make`.

## Usage
`gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]`

To check out all options and flags use `gfastats -h`.

## Assembly manipulation
**gfastats** allows extensive assembly manipulation at the sequence level. Manipulation is achieved using a set of *instructions* provided as an ordered list in a file to the option `-k` / `--swiss-army-knife`. The *instructions* are sequentially processed to generate the final output. Examples of *instructions* are:

```
JOIN contig1+ contig2+ 50 // creates a new gap of 50 bp between contig1 and contig2
SPLIT contig1+ contig2+ // splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them
REMOVE contig1 // removes contig1 from its scaffold, effectively splitting the scaffold in two
DELETE contig1:10-100 // deletes contig1 sequence between the coordinates provided (in bed format)
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
EXCISE contig1 50 // removes contig1 from the scaffold, leaving a 50 bp gap between the original sequences
INVERT contig1 // inverts contig1 sequence in place
RVCP contig1 // reverse-complement contig1 sequence in place
```

The *instructions* directly provide the list of edits that were introduced. The *instructions* could be from an automated tool or from manual annotation. See this wiki for a full list of *instructions*.

## Description
Please refer to **gfastats** paper for a complete description. Briefly, **gfastats** reads and stores any fasta<>fastq<>gfa[.gz] in gfa format. **gfastats** then builds a bidirected graph representation of the assembly using adjacency lists, where each node is a segment, and each edge is a gap (see figure below). The original sequence can be directly manipulated from the graph. Finally, walking the graph allows to generate different kinds of outputs, including manipulated assemblies and feature coordinates.

<p align="center">![alt gfastats assembly graph](images/graph.png)</p>

## How to cite
If you use **gfastats** in your work, please cite:

 
