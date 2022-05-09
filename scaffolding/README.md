### example data: bTaeGut2 Hifiasm (HiC) assembly
- hap1 contigs as GFA: `bTaeGut2.trim.HiC.hic.hap1.p_ctg.gfa` on [GenomeArk](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_hic_2.0/intermediates/hifiasm/)
- hap1 s1 AGP: `bTaeGut2_Saphyr_DLE1_3172351_bppAdjust_cmap_bTaeGut2_trim_HiC_hic_hap1_p_ctg_fasta_NGScontigs_HYBRID_SCAFFOLD.agp` on [GenomeArk](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_hic_2.0/intermediates/bionano_hap1/agp_fasta/)
- hap1 s1 edited AGP: 
- hap1 s2 AGP: `scaffolds_FINAL.original-coordinates.agp` on [GenomeArk](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Taeniopygia_guttata/bTaeGut2/assembly_vgp_hic_2.0/intermediates/salsa_hap1/bTaeGut2_hap1_s1.gfastats.rename_salsa/)

The starting files from hifiasm-HiC workflow are the hap1 & hap2 GFAs:

`bTaeGut2.hap1.gfa` and `bTaeGut2.hap2.gfa`

Convert GFA -> FASTA run bionano to obtain s1 AGPs. `bTaeGut2.hap1.fasta` into Bionano produces `bTaeGut2.hap1.s1.agp`, and same for hap2.

NOTE: IF Bionano is cutting, then fix the subseq lines. Bionano is not cutting in Galaxy, so do not need to run `sed` command on Galaxy assemblies.
````bash
# THIS IS NOT NEEDED FOR GALAXY ASSEMBLIES
cat bTaeGut2_hap1_s1.agp | sed 's/W\t\(.*\)_subseq_\([0-9]*\):\([0-9]*\)\t[0-9]*\t[0-9]*\t\(.\)/W\t\1\t\2\t\3\t\4/g' | sed 's/subseq_\([0-9]*\):\([0-9]*\)/subseq_\1_\2/g' > bTaeGut2_hap1_s1.edit.agp
````

##### UPDATE: MAY 3, 2022
Newer versions of gfastats append `_path` to path names, so the Bionano AGP must be processed accordingly. **This needs to happen even if Bionano is not cutting -- i.e. this needs to happen for Galaxy assemblies!**

an example of fixing the Bionano AGP to recognize `_path` in contig names:
````bash
awk '{OFS = "\t"}{if ($0 ~ /^#/) print $0 }{if ($6 ~ /h1*/) print $1,$2,$3,$4,$5,$6"_path",$7,$8,$9; if ($6 ~ /^[0-9]/) print $0}' bTaeGut2.hap1.s1.edit.agp > bTaeGut2.hap1.s1.edit.path.agp
````

Overlap s1 AGP onto c1/p1 GFA. `--discover` is so gfastats finds the paths in the GFA
````bash
gfastats bTaeGut2.hap1.gfa --discover -a bTaeGut2.hap1.s1.agp -o bTaeGut2.hap1.s1.gfa
````

Convert s1 GFA -> s1 FASTA, run salsa to obtain s2 AGP.
````bash
gfastats bTaeGut2.hap1.s1.gfa -o bTaeGut2.hap1.s1.gfastats.fasta
````
NOTE: IF Bionano is cutting, then subseq lines have colons in the names, so you need to remove those before SALSA
````bash
## Removing colons from bionano scaff names, because salsa doesn't like it
# THIS IS NOT NEEDED FOR GALAXY ASSEMBLIES
sed 's/:/_/g' bTaeGut2.hap1.s1.gfastats.fasta > bTaeGut2.hap1.s1.gfastats.nocolon.fasta
````

`bTaeGut2.hap1.s1.gfastats.fasta` into SALSA produces `bTaeGut2.hap1.s2.agp`

Overlap s2 AGP onto s1 GFA to create s2 GFA
````bash
cp <salsa_results_directory>/scaffolds_FINAL.original-coordinates.agp > ./bTaeGut2.hap1.s2.originalcoords.agp
gfastats bTaeGut2.hap1.s1.gfa -a bTaeGut2.hap1.s2.originalcoords.agp -o bTaeGut2.hap1.s2.gfa
````
If you want to convert this s2 GFA to s2 FASTA:
````bash
gfastats bTaeGut2.hap1.s2.gfa -o bTaeGut2.hap1.s2.gfastats.fasta
````
