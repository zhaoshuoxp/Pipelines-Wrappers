#!/bin/bash

# fisher-test-for-genomic-overlaps

#intersecting with bedtools
bedtools intersect -a $1 -b $3 > A_BG
bedtools intersect -a $2 -b $3 > B_BG

bedtools intersect -wa -a A_BG -b B_BG | uniq > A_B_BG
bedtools intersect -wa -a B_BG -b A_BG | uniq > B_A_BG
bedtools intersect -v -a A_BG -b A_B_BG | uniq > A_BG_noB
bedtools intersect -v -a B_BG -b B_A_BG | uniq > B_BG_noA
bedtools intersect -v -a $3 -b A_B_BG | uniq > BG_noA_noB

echo "Number of uniq A overlapping B in genomic background"
wc -l A_B_BG | cut -f1 -d ' '
echo "Number of uniq B overlapping A in genomic background"
wc -l B_A_BG | cut -f1 -d ' '
echo "Number of uniq A not overlapping B in genomic background"
wc -l A_BG_noB | cut -f1 -d ' '
echo "Number of uniq B not overlapping A in genomic background"
wc -l B_BG_noA | cut -f1 -d ' '
echo "Number of genomic background not overlapping A or B"

#writing and executing R script
echo "#!/usr/bin/env Rscript
A_B  <-
matrix(c($(wc -l A_B_BG | cut -f1 -d ' '), $(wc -l A_BG_noB | cut -f1 -d ' '), $(wc -l B_BG_noA | cut -f1 -d ' '), $(wc -l BG_noA_noB | cut -f1 -d ' ')),
       nrow = 2,
       dimnames = list(Guess = c(\"b+\", \"b-\"),
                       Truth = c(\"a+\", \"a-\")))
A_B.fisher <- fisher.test(A_B)
print(A_B.fisher)
str(A_B.fisher)" > script.r
chmod 775 script.r
./script.r

#removing intermediary files
rm A_BG
rm B_BG
rm A_B_BG
rm B_A_BG
rm A_BG_noB
rm B_BG_noA
rm BG_noA_noB
rm script.r
