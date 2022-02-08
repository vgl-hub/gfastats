# Instruction wiki

## JOIN

The JOIN instruction introduces a new gap of 50 bp between scaffold1 and scaffold2 with optional id gap1, effectively joining the two sequences into a new sequences named scaffold1 with an optional comment

```
JOIN contig1+ contig2+ 50 [gap1] [scaffold1] [this is a new scaffold]
```

## SPLIT

The SPLIT instruction splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them

```
SPLIT contig1+ contig2+
```

## REMOVE

The REMOVE instruction removes contig1, effectively splitting the scaffold in two if the contig is part of a scaffold

```
REMOVE contig1
```

## Yet to be implemented

```
DELETE contig1:10-100 // deletes contig1 sequence between the coordinates provided (in bed format)
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
EXCISE contig1 50 // removes contig1 from the scaffold, leaving a 50 bp gap between the original sequences
INVERT contig1 // inverts contig1 sequence in place
RVCP contig1 // reverse-complement contig1 sequence in place
```
