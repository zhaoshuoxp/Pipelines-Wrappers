#!/bin/bash
#####################################

# fisher-test-for-genes-overlaps
A=$1
B=$2
size=$(wc -l $A $B|tail -n 1 | awk '{print $1}')

# use whole genome gene set as background:
#wget https://dl.dropboxusercontent.com/s/9l3c5z9xsvdmofo/genelist.txt
# use size=A+B as background:
#sort -R genelist.txt |head -n $size > BG
#cat genelist.txt > BG
# use GWAS full catalog genes as background:
wget https://dl.dropboxusercontent.com/s/35zjtrvb1nj2vcy/gwas_gene.txt
cat gwas_gene.txt > BG
echo "background size:" 
wc -l BG

#intersecting with bedtools
cat $A BG |sort |uniq -d > A_BG
cat $B BG |sort |uniq -d > B_BG

cat A_BG B_BG |sort|uniq -d > A_B_BG
cat A_BG A_B_BG|sort|uniq -u > A_BG_noB
cat B_BG A_B_BG|sort|uniq -u > B_BG_noA
cat BG A_B_BG|sort|uniq -u > BG_noA_noB
cat BG A_BG|sort|uniq -u > BG_noA

echo "Number of uniq A overlapping B in background"
wc -l A_B_BG | awk '{print $1}'
echo "Number of uniq A not overlapping B in background"
wc -l A_BG_noB | awk '{print $1}'
echo "Number of uniq B not overlapping A in background"
wc -l B_BG_noA | awk '{print $1}'
echo "Number of background not overlapping A or B"
wc -l BG_noA_noB | awk '{print $1}'
echo "Number of background not overlapping A"
wc -l BG_noA | awk '{print $1}'

#writing and executing R script
echo "#!/usr/bin/env Rscript
A_B  <-
matrix(c($(wc -l A_B_BG | awk '{print $1}'), $(wc -l A_BG_noB | awk '{print $1}'), $(wc -l B_BG_noA | awk '{print $1}'), $(wc -l BG_noA_noB | awk '{print $1}')),
       nrow = 2,
       dimnames = list(Guess = c(\"b+\", \"b-\"),
                       Truth = c(\"a+\", \"a-\")))
A_B.fisher <- fisher.test(A_B)
print(A_B.fisher)
str(A_B.fisher)" > script.r
chmod 775 script.r
./script.r
rm script.r
#
## Hypergeometric-test-for-genes-overlaps
#echo "Hypergeometric p-value will be calculated with phyper in R as phyper(A-B overlap in BG, A no B in BG, total BG minus BG-A overlap, B no A in BG, lower.tail = FALSE, log.p = FALSE)"
##writing and executing R script
#echo "#!/usr/bin/env Rscript
#A_B.hypergeometric  <-
#phyper($(wc -l A_B_BG | awk '{print $1}'), $(wc -l A_BG_noB | awk '{print $1}'), $(wc -l BG_noA | awk '{print $1}'), $(wc -l B_BG_noA | awk '{print $1}'), lower.tail = FALSE, log.p = FALSE)
#A_B.hypergeometric
#str(A_B.hypergeometric)" > script.r
#
#chmod 775 script.r
#echo "Hypergeometric p-value:"
#./script.r
#rm script.r
#
#removing intermediary files
rm BG A_BG B_BG A_B_BG A_BG_noB B_BG_noA BG_noA_noB BG_noA
rm genelist.txt


################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################