testFiles/random2.gfa2 -o gfa
embedded
gfa gfa gfa
H	VN:Z:1.2
S	id2	TCAAGG
S	id3	CTTGATT
S	id1	ACCTT
S	id4	CATGACTC
S	id7	TGAATGAAA
L	id1	+	id2	-	3M
L	id2	-	id1	+	3M
J	id3	+	id4	-	5
J	id1	+	id2	+	3
J	id2	+	id3	-	2
J	id7	+	id7	+	5
P	id12	id1+;id2(1:3)+;id3-	*
P	path1	id1+	*
P	path2	id2+	*
P	path3	id3+	*
