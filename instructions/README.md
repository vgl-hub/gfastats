# Instruction wiki

Instructions are sequentially executed and each instruction is described by tab-separated columns.

## JOIN

The JOIN instruction introduces a new gap of 50 bp between `scaffold1` and `scaffold2` (two paths) with id `gap1`, effectively joining the two sequences into a new sequence with id `new_scaffold` and an optional comment.

```
JOIN    scaffold1+    scaffold2+    50  gap1    new_scaffold
JOIN    scaffold1(1:100)+    scaffold2(1:100)+    50  gap1    new_scaffold // optional subsetting
```

## SPLIT

The SPLIT instruction splits the scaffold containing `segment1` and `segment2`, effectively removing the existing gap between them. Two optional comments can be provided.

```
SPLIT   segment1+    segment2+    scaffold1   scaffold2 [this is a new scaffold1] [this is a new scaffold2]
```

## EXCISE

The EXCISE instruction removes segment1 from its scaffold, leaving it unplaced and adding a gap of the given size with id `gap1` between the original sequences

```
EXCISE  segment1  50  gap1
```

## REMOVE

The REMOVE instruction removes the paths involving the specified segment.

```
REMOVE  segment1
```

## EXCLUDE

The EXCLUDE instruction removes the specified path and all its components.

```
EXCLUDE  path1
```

## ERASE

The ERASE instruction trims off the sequence range specified from the given segment.

```
ERASE   segment1:10-100 // deletes segment1 sequence between the coordinates provided (in bed format)
```

## RVCP

The RVCP instruction reverse-complements segment1 sequence in place

```
RVCP    segment1
```

## INVERT

The INVERT instruction inverts segment1 sequence in place

```
INVERT  segment1
```

## Yet to be implemented

```
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
```
