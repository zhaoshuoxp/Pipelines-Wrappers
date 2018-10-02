#!/bin/bash

#get gwas plus CARDIoGRAMplusC4D SNPs in BED format
#wget http://www.genome.gov/admin/gwascatalog.txt
awk -v OFS="\t" -F"\t" '{print $12,$13,$13+1,$22,$8}' /home/quanyi/SNP_dataset/gwascatalog.txt | awk -F"\t" '{if ($3!=1){print $0}}' > GwasCatalog.bed
awk -v OFS="\t" '{print $1,$2,$3,$4,"CardiogramPlusC4D"}' /home/quanyi/SNP_dataset/CARDIoGRAMplusC4D/CARDIOGRAMplusC4DleadSNPs.bed > CARDIOGRAMC4Dplusnovel.tmp
sed -i 's/^chr//g' CARDIOGRAMC4Dplusnovel.tmp
cat CARDIOGRAMC4Dplusnovel.tmp >> GwasCatalog.bed
rm CARDIOGRAMC4Dplusnovel.tmp

#writing main to perform grep on gwascatalog.txt.cut
echo "
#grep categories
for var in \"\${@:1:\$#-1}\"
do
	echo Received: \$var
grep -i  \"\$var\" GwasCatalog.bed > \"\$var\".gwascatalog.bed 
echo Done: \$var
done

#print stats
echo "Gwas Catalog number of SNP-phenotype associations:"
wc -l GwasCatalog.bed
echo "Gwas Catalog number of SNP-phenotype associations per category:"

for var in \"\${@:1:\$#-1}\"
do
	echo Phenotype: \$var
wc -l \"\$var\".gwascatalog.bed 
done

#bed files - cut, sort, and uniq
for var in \"\${@:1:\$#-1}\"
do
cut -f1-3 \"\$var\".gwascatalog.bed > \"\$var\".gwascatalog.bed.cut
sort -k1,1V -k2,2n \"\$var\".gwascatalog.bed.cut > \"\$var\".gwascatalog.bed.cut.sort
uniq \"\$var\".gwascatalog.bed.cut.sort > \"\$var\".gwascatalog.bed.cut.sort.uniq 
rm \"\$var\".gwascatalog.bed.cut
rm \"\$var\".gwascatalog.bed.cut.sort
done

#new stats
echo "Gwas Catalog number of SNP-phenotype associations per category AFTER REMOVING DUPLICATES:"

for var in \"\${@:1:\$#-1}\"
do
	echo Phenotype: \$var
wc -l \"\$var\".gwascatalog.bed.cut.sort.uniq

done

#substitute 23 24 with X Y, add chr
for var in \"\${@:1:\$#-1}\"
do
	echo Converting Phenotype: \$var
sed -i 's/^23/X/g' \"\$var\".gwascatalog.bed.cut.sort.uniq
sed -i 's/^24/Y/g' \"\$var\".gwascatalog.bed.cut.sort.uniq
sed  's/^/chr/g' \"\$var\".gwascatalog.bed.cut.sort.uniq >  \"\$var\".gwascatalog.bed.cut.sort.uniq.chrXY
rm \"\$var\".gwascatalog.bed.cut.sort.uniq
done

" > main.sh

#calling main
chmod 775 main.sh
./main.sh "$@"

#removing files
rm main.sh

for last; do true; done
echo $last

echo Input phenotypes:
echo "${@:1:$#-1}"

#overlapping  
for var in "${@:1:$#-1}"
do
	echo Overlapping Phenotype SNPs with input bed: $var
bedtools intersect -a "$var".gwascatalog.bed.cut.sort.uniq.chrXY -b $last > "$var".gwascatalog.bed.cut.sort.uniq.overlap
done

for var in "${@:1:$#-1}"
do
	echo Number of Overlapping Phenotype SNPs with input bed: $var
wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap
done

for var in "${@:1:$#-1}"
do
	echo Number of Overlapping Phenotype SNPs with input bed: $var
wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap | cut -f1 -d ' '
done
#get the size of hg19
hg19=$(cat $len_hg19 | awk '{ sum+=($2)} END {print sum}')
echo Human Genome size version hg19: $hg19

#bedtools merge on input bed
sort -k1,1V -k2,2n $last > tmp
bedtools merge -i tmp -c 1 -o count > $last
wc -l $last
rm tmp
 
#calculate coverage
cov=$(cat $last | awk '{ sum+=($3-$2)} END {print sum}')
echo Coverage of BED file $cov
fra=$(cat $last | awk '{ sum+=($3-$2)} END {print sum/"'"$hg19"'"}')
echo Fraction of hg19 $fra

#make overlaps with input bed ranges
for var in "${@:1:$#-1}"
do
	echo Overlapping Phenotype SNPs with input bed to calculate coverage: $var
bedtools intersect -wb -a "$var".gwascatalog.bed.cut.sort.uniq.chrXY -b $last > "$var".gwascatalog.bed.cut.sort.uniq.overlap.input.int
cut -f4-6 "$var".gwascatalog.bed.cut.sort.uniq.overlap.input.int > "$var".gwascatalog.bed.cut.sort.uniq.overlap.input.int.cut
done

# prepare for fold change cal 
selected=$(wc -l *uniq*cut|tail -n 1|cut -f 3 -d ' ')
selected_total=$(wc -l *uniq*chrXY|tail -n 1|cut -f 3 -d ' ')
sed -i 1d GwasCatalog.bed
awk -v OFS="\t" '{print "chr"$1,$2,$3}' GwasCatalog.bed > GwasCatalog.bed1
all=$(wc -l GwasCatalog.bed1|cut -f 1 -d" ")
sub_total=$(intersectBed -a GwasCatalog.bed1 -b "$last" |wc -l|cut -d" " -f1)

#create  R script
touch script.R
echo "#!/usr/bin/env Rscript
existingDF <- as.data.frame(matrix(seq(4),nrow=1,ncol=4))" >script.R
#calculate coverage and fraction per category, load into R script
i=1
for var in "${@:1:$#-1}"
do
let i=i+1
#removing input spaces and '
var2a=$(echo $var | tr -d ' ')
var2=${var2a//[^[:alnum:]]/}
echo Input terms no spaces: $var2

var3="$(cat "$var".gwascatalog.bed.cut.sort.uniq.overlap.input.int.cut | awk '{ sum+=($3-$2)} END {print sum}')"
echo Coverage of the input bed file that overlaps GWAS category: $var3

fra=$(cat "$var".gwascatalog.bed.cut.sort.uniq.overlap.input.int.cut | awk '{ sum+=($3-$2)} END {print sum/"'"$hg19"'"}')
echo Fraction of hg19 $fra

echo "print(\"$var\")" >> script.R
echo  "dbinom ($(wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap | cut -f1 -d ' '), $(wc -l "$var".gwascatalog.bed | cut -f1 -d ' '), "$fra")" >> script.R
#calculating fold change

echo Fold change: $var2
overlap="$(wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap | cut -f1 -d ' ')"
total="$(wc -l "$var".gwascatalog.bed | cut -f1 -d ' ')"
echo overlap: $overlap
echo total: $total
echo fold:
# fold change to all catalogs
fold="$(awk 'BEGIN {print (("'"$overlap"'"/"'"$sub_total"'")/("'"$total"'"/"'"$all"'"))}')"
# fold change to selected catalogs
#fold="$(awk 'BEGIN {print (("'"$overlap"'"/"'"$selected"'")/("'"$total"'"/"'"$selected_total"'"))}')"
echo $fold

echo  "x<-dbinom ($(wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap | cut -f1 -d ' '), $(wc -l "$var".gwascatalog.bed | cut -f1 -d ' '), "$fra")
y<-"$fold"
name<-\""$var2"\"
s<-$(wc -l "$var".gwascatalog.bed | cut -f1 -d ' ')
z<-c(name,x,y,s)
existingDF <- rbind(existingDF,z)
existingDF">>script.R
done

#finishing R script

echo "existingDF<-existingDF[-1,]
existingDF[,2:4]<-sapply(existingDF[,2:4], as.numeric)
sapply(existingDF, mode)
existingDF
existingDF<-transform(existingDF, V2=-log(V2))
existingDF

data<-existingDF
data
rownames(data) <- data\$V1
data<-data[,2:4]

#adding categories
k<-dim (data)
rep <-rep(\"Other\", k[1])
data\$V5 <- rep

data[rownames(data) == \"Schizophrenia\",]\$V5<-\"Brain\"
data[rownames(data) == \"Bipolardisorder\",]\$V5<-\"Brain\"
data[rownames(data) == \"Atherosclerosis\",]\$V5<-\"Cardiovascular\"
data[rownames(data) == \"Ulcerativecolitis\",]\$V5<-\"Chronic Inflammatory\"
data[rownames(data) == \"CardiogramplusC4D\",]\$V5<-\"Cardiovascular\"
data[rownames(data) == \"CoronaryHeart\",]\$V5<-\"Cardiovascular\"
data[rownames(data) == \"Myocardialinfarction\",]\$V5<-\"Cardiovascular\"
data[rownames(data) == \"Multiplesclerosis\",]\$V5<-\"Brain\"
data[rownames(data) == \"Parkinsonsdisease\",]\$V5<-\"Brain\"
data[rownames(data) == \"Alzheimersdisease\",]\$V5<-\"Brain\"
data[rownames(data) == \"Lupus\",]\$V5<-\"Chronic Inflammatory\"
data[rownames(data) == \"Prostatecancer\",]\$V5<-\"Cancer\"
data[rownames(data) == \"Pancreaticcancer\",]\$V5<-\"Cancer\"
data[rownames(data) == \"Breastcancer\",]\$V5<-\"Cancer\"
data[rownames(data) == \"CoronaryArtery\",]\$V5<-\"Cardiovascular\"
data[rownames(data) == \"Coronaryarterycalcification\",]\$V5<-\"Cardiovascular\"

#removing 0s
data[, 1:3] <- sapply(data[,1:3], as.numeric)
row_sub = apply(data[,1:3], 1, function(y) all(y != 0))
row_sub
data<-data[row_sub,]
data
data[, 1:3] <- sapply(data[,1:3], as.numeric)
colnames(data)<-c(\"LogP\", \"FC\", \"Phenotype SNPs\", \"Category\")
data

#making ggplot2 graph
library(ggplot2)
library(wesanderson)
library(directlabels)
ymax<-max(data\$LogP,na.rm = TRUE)
ymin<-min(data\$LogP,na.rm = TRUE)
xmax<-max(data\$FC,na.rm = TRUE)
xmin<-min(data\$FC,na.rm = TRUE)

p<- ggplot(data, aes(x=data\$FC, y=data\$LogP,label=row.names(data))) + geom_point(shape=19, alpha=1/8, color=\"red\", aes(size=data\$\"Phenotype SNPs\"), max_size=max(data\$\"Phenotype SNPs\")) + xlab(\"Fold change\") + ylab(\"-log P-value\") + ggtitle (\"GWAS SNPs enrichment - binomial test\") + geom_dl(aes(label=row.names(data)), method=list(\"first.bumpup\"), col=\"blue\", alpha=1/2)+xlim(xmin-1, xmax)
pdf(\"output.pdf\",width=10, height=8)
print(p+ geom_dl(aes(colour = data\$\"Category\"), method=list(\"first.bumpup\")) + scale_colour_hue(name=\"Category\") + labs(size=\"Phenotype SNPs\", color=\"Category\") + scale_size(range = c(0,50)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face=\"bold\")) + scale_color_manual(values = c(wes_palette(\"Cavalcanti1\"), wes_palette(\"Royal1\"), wes_palette(\"GrandBudapest2\"), wes_palette(\"Royal2\"), wes_palette(\"Darjeeling2\"), wes_palette(\"Zissou1\")))) 
dev.off()">>script.R


chmod 775 script.R
./script.R
rm script.R
rm *.gwascatalog* GwasCatalog.bed GwasCatalog.bed1
