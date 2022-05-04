The starting files from hifiasm-HiC workflow are the hap1 & hap2 GFAs:

`bColStr4.hap1.gfa` and `bColStr4.hap2.gfa`

Convert GFA -> FASTA run bionano to obtain s1 AGPs. `bColStr4.hap1.fasta` into Bionano produces `bColStr4.hap1.s1.agp`, and same for hap2.

NOTE: IF Bionano is cutting, then fix the subseq lines. Bionano is not cutting in Galaxy, so do not need to run `sed` command on Galaxy assemblies.
````bash
# THIS IS NOT NEEDED FOR GALAXY ASSEMBLIES
sed 's/W\t\(.*\)_subseq_\([0-9]*\):\([0-9]*\)\t[0-9]*\t[0-9]*\t\(.\)/W\t\1\t\2\t\3\t\4/g' bTaeGut2_hap1_s1.agp > bTaeGut2_hap1_s1.edit.agp
````

##### UPDATE: MAY 3, 2022
Newer versions of gfastats append `_path` to path names, so the Bionano AGP must be processed accordingly. **This needs to happen even if Bionano is not cutting -- i.e. this needs to happen for Galaxy assemblies!**

an example of fixing the Bionano AGP to recognize `_path` in contig names:
````bash
sed 's/\(h*tg[0-9]*.\)/\1_path/g' bTaeGut2.hap1.s1.edit.agp | sed 's/_path_obj/_obj/g' | sed 's/_path_subseq_/_subseq_/g' | sed 's/:/_/g' > bTaeGut2.hap1.s1.edit.path.agp
````

Overlap s1 AGP onto c1/p1 GFA. `--discover` is so gfastats finds the paths in the GFA
````bash
gfastats bColStr4.hap1.gfa --discover -a bColStr4.hap1.s1.agp -o bColStr4.hap1.s1.gfa
````

Convert s1 GFA -> s1 FASTA, run salsa to obtain s2 AGP.
````bash
gfastats bColStr4.hap1.s1.gfa -o bColStr4.hap1.s1.gfastats.fasta
````
NOTE: IF Bionano is cutting, then subseq lines have colons in the names, so you need to remove those before SALSA
````bash
## Removing colons from bionano scaff names, because salsa doesn't like it
# THIS IS NOT NEEDED FOR GALAXY ASSEMBLIES
sed 's/:/_/g' bColStr4.hap1.s1.gfastats.fasta > bColStr4.hap1.s1.gfastats.nocolon.fasta
````

`bColStr4.hap1.s1.gfastats.fasta` into SALSA produces `bColStr4.hap1.s2.agp`

Overlap s2 AGP onto s1 GFA to create s2 GFA
````bash
cp <salsa_results_directory>/scaffolds_FINAL.original-coordinates.agp > ./bColStr4.hap1.s2.originalcoords.agp
gfastats bColStr4.hap1.s1.gfa -a bColStr4.hap1.s2.originalcoords.agp -o bColStr4.hap1.s2.gfa
````
If you want to convert this s2 GFA to s2 FASTA:
````bash
gfastats bColStr4.hap1.s2.gfa -o bColStr4.hap1.s2.gfastats.fasta
````
