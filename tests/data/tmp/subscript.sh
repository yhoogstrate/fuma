

for chr in $( cut -f 1 test_ParseBED.TestParseBED.bed | sort | uniq )
do
	echo ">> "$chr" <<"
	grep $chr "test_ParseBED.TestParseBED.bed" | head -n 2 >> bla.txt
done