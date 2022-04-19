The starting files from hifiasm-HiC workflow are the hap1 & hap2 GFAs:

`bColStr4.hap1.gfa` and `bColStr4.hap2.gfa`

Convert GFA -> FASTA run bionano to obtain s1 AGPs. `bColStr4.hap1.fasta` into Bionano produces `bColStr4.hap1.s1.agp`, and same for hap2.

NOTE: IF Bionano is cutting, then fix the subseq lines. Bionano is not cutting in Galaxy, so do not need to run `sed` command on Galaxy assemblies.
````bash
## THIS IS NOT NEEDED FOR GALAXY ASSEMBLIES
sed 's/W\t\(.*\)_subseq_\([0-9]*\):\([0-9]*\)\t[0-9]*\t[0-9]*\t\(.\)/W\t\1\t\2\t\3\t\4/g' bTaeGut2_hap1_s1.agp > bTaeGut2_hap1_s1.edit.agp
````

Overlap s1 AGP onto c1/p1 GFA. `--discover` is so gfastats finds the paths in the GFA
````bash
$GFASTATS bColStr4.hap1.gfa --discover -a bColStr4.hap1.s1.agp -o bColStr4.hap1.s1.gfa
````

Convert s1 GFA -> s1 FASTA, run salsa to obtain s2 AGP. `--line-length 60` is telling gfastats to output FASTA wrapped at 60th char
````bash
$GFASTATS bColStr4.hap1.s1.gfa --line-length 60 -o bColStr4.hap1.s1.gfastats.fasta
# Need to remove colons from bionano scaff names, because salsa doesn't like it
sed 's/:/_/g' bColStr4.hap1.s1.gfastats.fasta > bColStr4.hap1.s1.gfastats.nocolon.fasta

<run salsa on bColStr4.hap1.s1.gfastats.nocolon.fasta to obtain bColStr4.hap1.s2.agp>
````

Overlap s2 AGP onto s1 GFA, convert s2 GFA to s2 FASTA
````bash
cp <salsa_results_directory>/scaffolds_FINAL.original-coordinates.agp > ./bColStr4.hap1.s2.originalcoords.agp
$GFASTATS bColStr4.hap1.s1.gfa -a bColStr4.hap1.s2.originalcoords.agp -o bColStr4.hap1.s2.gfa
$GFASTATS bColStr4.hap1.s2.gfa --line-length 60 -o bColStr4.hap1.s2.gfastats.fasta
````
