# gfastats
A single fast and exhaustive tool to compute summary statistics and simultaneously manipulate \*fast\* (fasta, fastq, gfa [.gz]) files.
It also allows seamless fasta<>fastq<>gfa[.gz] conversion.

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

Coordinates ans sizes of scaffolds, contigs and gaps can be outputted with the option `--out-coord` and `--out-size`.

Additionally, input can be filtered using scaffold lists or `bed` coordinate files with the options `--include-bed` and `--exclude-bed`.

Importantly, the filtered input can be outputted in any \*fast\* (fasta, fastq, gfa [.gz]) format.

## Installation
Either download one of the releases or `git clone https://github.com/vgl-hub/gfastats.git` and `make`.

## Usage
`gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]`
To check out all option use `gfastats -h`.

##Description
 
