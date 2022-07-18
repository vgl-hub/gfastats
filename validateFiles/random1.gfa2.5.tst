testFiles/random1.gfa2 -o gfa2
embedded
H	VN:Z:2.0
S	Header1.1	5	CGacT
S	Header2.1	3	CGA
S	Header2.3	1	T
S	Header3.1	3	TGA
S	Header3.3	2	AT
S	Header3.5	2	CT
S	Header4.2	12	TTCCTcgCACtC
S	Header5.1	12	AACTCGATCACG
G	Header2.2	Header2.1+	Header2.3+	1
G	Header3.2	Header3.1+	Header3.3+	1
G	Header3.4	Header3.3+	Header3.5+	1
G	Header3.6	Header3.5+	Header3.5-	1
G	Header4.1	Header4.2+	Header4.2+	3
G	Header5.2	Header5.1+	Header5.1-	3
O	Header1	Header1.1+
O	Header2	Header2.1+ Header2.2 Header2.3+
O	Header3	Header3.1+ Header3.2 Header3.3+ Header3.4 Header3.5+ Header3.6
O	Header4	Header4.1 Header4.2+
O	Header5	Header5.1+ Header5.2
