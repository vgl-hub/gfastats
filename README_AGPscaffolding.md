the starting files from hifiasm-HiC workflow:
````bash
lrwxrwxrwx  1 labueg jarv   54 Apr 11 15:17 bTaeGut2.trim.HiC.hic.hap1.p_ctg.gfa -> ../2d_hifiasm_hic/bTaeGut2.trim.HiC.hic.hap1.p_ctg.gfa
lrwxrwxrwx  1 labueg jarv   54 Apr 11 15:17 bTaeGut2.trim.HiC.hic.hap2.p_ctg.gfa -> ../2d_hifiasm_hic/bTaeGut2.trim.HiC.hic.hap2.p_ctg.gfa
````
Convert GFA -> FASTA, run bionano to obtain s1 AGPs (here I am just using the agp files from the hifiasm-HiC scaffolding run, as there was no purging. so c1 hap1/2 -> s1 hap1/2)

Copy (or link) the s1 AGPs and fix the subseq lines
````bash
cp ../3a_bionano_hap1/agp_fasta/bTaeGut2_Saphyr_DLE1_3172351_bppAdjust_cmap_bTaeGut2_trim_HiC_hic_hap1_p_ctg_fasta_NGScontigs_HYBRID_SCAFFOLD.agp ./bTaeGut2_hap1_s1.agp
sed 's/W\t\(.*\)_subseq_\([0-9]*\):\([0-9]*\)\t[0-9]*\t[0-9]*\t\(.\)/W\t\1\t\2\t\3\t\4/g' bTaeGut2_hap1_s1.agp > bTaeGut2_hap1_s1.edit.agp
````

Overlap s1 AGP onto c1/p1 GFA. `--discover` is so gfastats finds the paths in the GFA
````bash
$GFASTATS bTaeGut2.trim.HiC.hic.hap1.p_ctg.gfa --discover -a bTaeGut2_hap1_s1.edit.agp -o bTaeGut2.hap1.s1.gfa
````

Convert s1 GFA -> s1 FASTA, run salsa to obtain s2 AGP
````bash
$GFASTATS bTaeGut2.hap1.s1.gfa --line-length 60 -o bTaeGut2.hap1.s1.gfastats.fasta
# optional: change bionano scaff names to match w normal pipeline names
# sed 's/-/_/g' bTaeGut2.hap1.s1.gfastats.fasta > bTaeGut2.hap1.s1.gfastats.renamed.fasta
# NOT optional: remove colons from bionano scaff names, because salsa doesn't like it
sed 's/:/_/g' bTaeGut2.hap1.s1.gfastats.fasta > bTaeGut2.hap1.s1.gfastats.nocolon.fasta
sh $STORE/programs/SALSA_scripts/_submit_salsa_2.2.sh ../bTaeGut2.hap1.s1.gfastats.fasta $BTAEGUT2/genomic_data/arima/ vgl 32
````

Overlap s2 AGP onto s1 GFA, convert s2 GFA to s2 FASTA
````bash
cp [salsa_results_directory]/scaffolds_FINAL.original-coordinates.agp > ./bTaeGut2_hap1_s2_originalcoords.agp
$GFASTATS bTaeGut2.hap1.s1.gfa -a bTaeGut2_hap1_s2_originalcoords.agp -o bTaeGut2.hap1.s2.gfa
$GFASTATS bTaeGut2.hap1.s2.gfa --line-length 60 -o bTaeGut2.hap1.s2.gfastats.fasta
````
