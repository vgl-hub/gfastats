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

The EXCISE instruction removes segment1 from its scaffold, leaving it unplaced and adding a gap of 50bp with id `gap1` between the original sequences

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

The RVCP instruction reverse-complements path1 or segment1 sequence in place

```
RVCP    path1/segment1
```

## INVERT

The INVERT instruction inverts segment1 sequence in place

```
INVERT  segment1
```

## RESIZE

The RESIZE instruction resizes the size of gap1 to 50 bp

```
RESIZE  gap1    50
```

## MASK

The MASK instruction masks with 50 Ns a portion of a path, effectively adding a gap in the corresponding segment of optional size 5. If size is not provided, the masked size is used

```
MASK  path1 10  60  [5]
```

## CLEAVE

The CLEAVE instruction breaks the specified segment at the given position generating segment2 and segment3, optionally connected by an edge

```
CLEAVE  segment1 50  segment2 segment3 [edge1]
```

## Yet to be implemented

```
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
```
