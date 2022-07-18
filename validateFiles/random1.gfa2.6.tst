testFiles/random1.gfa2 -o gfa
embedded
gfa gfa gfa
H	VN:Z:1.2
S	Header1.1	CGacT
S	Header2.1	CGA
S	Header2.3	T
S	Header3.1	TGA
S	Header3.3	AT
S	Header3.5	CT
S	Header4.2	TTCCTcgCACtC
S	Header5.1	AACTCGATCACG
J	Header2.1	+	Header2.3	+	1
J	Header3.1	+	Header3.3	+	1
J	Header3.3	+	Header3.5	+	1
J	Header3.5	+	Header3.5	-	1
J	Header4.2	+	Header4.2	+	3
J	Header5.1	+	Header5.1	-	3
P	Header1	Header1.1+	*
P	Header2	Header2.1+;Header2.3+	*
P	Header3	Header3.1+;Header3.3+;Header3.5+;	*
P	Header4	;Header4.2+	*
P	Header5	Header5.1+;	*
