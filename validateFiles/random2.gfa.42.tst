testFiles/random2.gfa -o gfa2
embedded
H	VN:Z:2.0
S	11	5	ACCTT	LN:i:5	QL:Z:?@97?
S	12	6	TCAAGG	LN:i:6	QL:Z:@6?84@
S	13	7	CTTgaTT	LN:i:7	QL:Z:>=?@877
E	edge0	11	+	12	-	4M
E	edge1	12	-	13	+	5M
E	edge2	11	+	13	+	3M
G	gap0	11+	13-	5	SC:i:1
G	gap1	13-	12+	3	SC:i:1
O	14	11+ gap0 13- gap1 12+
O	15	11+ edge0 12- edge1 13+
