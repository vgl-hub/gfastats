# Instruction wiki

##JOIN

```
JOIN contig1+ contig2+ 50 [scaffold1] [this is a new scaffold] // introduces a new gap of 50 bp between scaffold1 and scaffold2, effectively joining the two sequences into a new sequences named scaffold1 with an optional comment
```

##SPLIT

```
SPLIT contig1+ contig2+ // splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them
```

##Yet to be implemented

```
REMOVE contig1 // removes contig1 from its scaffold, effectively splitting the scaffold in two
DELETE contig1:10-100 // deletes contig1 sequence between the coordinates provided (in bed format)
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
EXCISE contig1 50 // removes contig1 from the scaffold, leaving a 50 bp gap between the original sequences
INVERT contig1 // inverts contig1 sequence in place
RVCP contig1 // reverse-complement contig1 sequence in place
```
