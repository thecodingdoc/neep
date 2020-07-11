#!/bin/sh
PVAL=0.01

if [ ! -e ../neep ]
then
    echo "Build neep first. In distribution directory, run 'make'."
    exit 1
fi

if [ ! -e test_output.txt ]
then
    ../neep -c test_clinical.csv -e test_expression.csv -o test_output.txt -n 20000 -t 0.15
else
    echo **note: using previously generated 'test_output.txt'.
    echo **note: delete if a re-run is desired.
fi

perl -F"\t" -lane "print join(\" \",\$F[0],\$F[5]) if \$F[3]<$PVAL" correct_test_output.txt | sort -k 1 | cut -d" " -f1,2 > correct.txt
perl -F"\t" -lane "print join(\" \",\$F[0],\$F[5]) if \$F[3]<$PVAL" test_output.txt | sort -k 1 | cut -d" " -f1,2 > test.txt

if join correct.txt test.txt | perl -lane 'BEGIN{ $good = 1;}' -e'if ($F[1] ne $F[2]) { $good=0; print; }' -e 'END{ $good ? exit 0 : exit 1;}'
then
    N=$(join correct.txt test.txt | wc -l)
    echo "success ($N calls examined at p value $PVAL)";
else
    echo "FAILURE - test output materially different from correct_test_output.txt"
fi
