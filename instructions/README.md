# Instruction wiki

## JOIN

The JOIN instruction introduces a new gap of 50 bp between scaffold1 and scaffold2 with id gap1, effectively joining the two sequences into a new sequences named scaffold1 with an optional comment.

```
JOIN    contig1+    contig2+    50  gap1    scaffold1   [this is a new scaffold]
```

## SPLIT

The SPLIT instruction splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them. Two optional comments can be provided.

```
SPLIT   contig1+    contig2+    scaffold1   scaffold2 [this is a new scaffold1] [this is a new scaffold2]
```

## EXCISE

The EXCISE instruction removes contig1 from its scaffold, leading leaving it unplaced and adding a gap of the given size with optional id gap1 between the original sequences

```
EXCISE  contig1  50  [gap1] // new 50 bp gap
```

## REMOVE

The REMOVE instruction removes contig1 from the segment set. If it is part of a path, the path is also removed.

```
REMOVE  contig1
```

## ERASE

The ERASE instruction trims off the sequence range specified from the given segment.

```
ERASE   contig1:10-100 // deletes contig1 sequence between the coordinates provided (in bed format)
```

## RVCP

The RVCP instruction reverse-complements contig1 sequence in place

```
RVCP    contig1
```

## INVERT

The INVERT instruction inverts contig1 sequence in place

```
INVERT  contig1
```

## Yet to be implemented

```
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
```
