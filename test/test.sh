#!/bin/sh

neep -c test_clinical.csv -e test_expression.csv -o test_output.txt -n 20000 -t 0.15

echo "complete";
