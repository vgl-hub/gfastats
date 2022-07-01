testFiles/random2.gfa2 -o gfa2
embedded
H	VN:Z:2.0
S	id2	6	TCAAGG
S	id3	7	CTTGATT
S	id1	5	ACCTT
S	id4	8	CATGACTC
S	id7	9	TGAATGAAA
E	id1	+	id2	-	3M
E	id2	-	id1	+	3M
G	id5	id3+	id4-	5
G	id6	id1+	id2+	3
G	id8	id2+	id3-	2
G	id9	id7+	id7+	5
O	id12	id1+ id6 id2(1:3)+ id8 id3-
O	path1	id1+
O	path2	id2+
O	path3	id3+
