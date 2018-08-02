#!/bin/bash
#download gwas catalog and create bed file with "chr position position+1 proxy_gene phenotype"
#echo "
##download gwas catalog 
#awk -F\"\t\" '{if (\$12!=\"\") print \$12\"\t\"\$13\"\t\"\$15\"\t\"\$8}' /home/quanyi/SNP_dataset/gwascatalog.txt > tmp
#awk -F\"\t\" '{print \$1\"\t\"\$2\"\t\"\$2+1\"\t\"\$3\"\t\"\$4}' tmp > GwasCatalog.bed
#rm tmp" > GwasCatalog2Bed.sh
#chmod 775 GwasCatalog2Bed.sh
#./GwasCatalog2Bed.sh
#rm GwasCatalog2Bed.sh

#adding CardiogramPlusC4D to GWASCatalog
#awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t","CardiogramPlusC4D"}' /home/quanyi/SNP_dataset/CARDIoGRAMplusC4D/CARDIOGRAMplusC4DleadSNPs.bed  > CARDIOGRAMC4Dplusnovel.txt.tmp
#sed -i 's/^chr//g' CARDIOGRAMC4Dplusnovel.txt.tmp
#cat CARDIOGRAMC4Dplusnovel.txt.tmp >> GwasCatalog.bed
#rm CARDIOGRAMC4Dplusnovel.txt.tmp

# bugs above, some pos=NA appear, use below one
awk -v OFS="\t" -F"\t" '{print $12,$13,$13+1,$15,$8}' /home/quanyi/SNP_dataset/gwascatalog.txt | awk -F"\t" '{if ($3!=1){print $0}}' > GwasCatalog.bed
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
echo "#!/usr/bin/Rscript" > script.R
echo "existingDF <- as.data.frame(matrix(seq(4),nrow=1,ncol=4))" >>script.R
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

echo  "x<-dbinom ($(wc -l "$var".gwascatalog.bed.cut.sort.uniq.overlap | cut -f1 -d ' '), $(wc -l "$var".gwascatalog.bed | cut -f1 -d ' '), "$fra")" >> script.R
echo "y<-"$fold"">> script.R
echo "name<-\""$var2"\"">>script.R
echo "s<-$(wc -l "$var".gwascatalog.bed | cut -f1 -d ' ')">>script.R
echo "z<-c(name,x,y,s)">>script.R
echo "existingDF <- rbind(existingDF,z)">>script.R
echo "existingDF">>script.R
done

#finishing R script

echo "existingDF<-existingDF[-1,]">>script.R
echo "existingDF[,2:4]<-sapply(existingDF[,2:4], as.numeric)">>script.R
echo "sapply(existingDF, mode)">>script.R
echo "existingDF">>script.R
echo "existingDF<-transform(existingDF, V2=-log(V2))" >> script.R
echo "existingDF">>script.R

echo "data<-existingDF" >>script.R
echo "data">>script.R
echo "rownames(data) <- data\$V1" >>script.R
echo "data<-data[,2:4]">>script.R

#adding categories
echo "k<-dim (data)" >>script.R
echo "rep <-rep(\"Other\", k[1])">>script.R
echo "data\$V5 <- rep">>script.R

echo "data[rownames(data) == \"Schizophrenia\",]\$V5<-\"Brain\"">>script.R
echo "data[rownames(data) == \"Bipolardisorder\",]\$V5<-\"Brain\"">>script.R
echo "data[rownames(data) == \"Atherosclerosis\",]\$V5<-\"Cardiovascular\"">>script.R
echo "data[rownames(data) == \"Ulcerativecolitis\",]\$V5<-\"Chronic Inflammatory\"">>script.R
echo "data[rownames(data) == \"CardiogramplusC4D\",]\$V5<-\"Cardiovascular\"">>script.R
echo "data[rownames(data) == \"CoronaryHeart\",]\$V5<-\"Cardiovascular\"">>script.R
echo "data[rownames(data) == \"Myocardialinfarction\",]\$V5<-\"Cardiovascular\"">>script.R
echo "data[rownames(data) == \"Multiplesclerosis\",]\$V5<-\"Brain\"">>script.R
echo "data[rownames(data) == \"Parkinsonsdisease\",]\$V5<-\"Brain\"">>script.R
echo "data[rownames(data) == \"Alzheimersdisease\",]\$V5<-\"Brain\"">>script.R
echo "data[rownames(data) == \"Lupus\",]\$V5<-\"Chronic Inflammatory\"">>script.R
echo "data[rownames(data) == \"Prostatecancer\",]\$V5<-\"Cancer\"">>script.R
echo "data[rownames(data) == \"Pancreaticcancer\",]\$V5<-\"Cancer\"">>script.R
echo "data[rownames(data) == \"Breastcancer\",]\$V5<-\"Cancer\"">>script.R
echo "data[rownames(data) == \"CoronaryArtery\",]\$V5<-\"Cardiovascular\"">>script.R
echo "data[rownames(data) == \"Coronaryarterycalcification\",]\$V5<-\"Cardiovascular\"">>script.R

#removing 0s
echo "data[, 1:3] <- sapply(data[,1:3], as.numeric)">>script.R
echo "row_sub = apply(data[,1:3], 1, function(y) all(y != 0))">>script.R
echo "row_sub">>script.R
echo "data<-data[row_sub,]">>script.R
echo "data">>script.R
echo "data[, 1:3] <- sapply(data[,1:3], as.numeric)">>script.R
echo "colnames(data)<-c(\"LogP\", \"FC\", \"Phenotype SNPs\", \"Category\")">>script.R
echo "data">>script.R

#making ggplot2 graph
echo "library(ggplot2)" >> script.R
echo "library(wesanderson)">>script.R
echo "library(directlabels)">>script.R
echo "ymax<-max(data\$LogP,na.rm = TRUE)">>script.R
echo "ymin<-min(data\$LogP,na.rm = TRUE)">>script.R
echo "xmax<-max(data\$FC,na.rm = TRUE)">>script.R
echo "xmin<-min(data\$FC,na.rm = TRUE)">>script.R

echo "p<- ggplot(data, aes(x=data\$FC, y=data\$LogP,label=row.names(data))) + geom_point(shape=19, alpha=1/8, color=\"red\", aes(size=data\$\"Phenotype SNPs\"), max_size=max(data\$\"Phenotype SNPs\")) + xlab(\"Fold change\") + ylab(\"-log P-value\") + ggtitle (\"GWAS SNPs enrichment - binomial test\") + geom_dl(aes(label=row.names(data)), method=list(\"first.bumpup\"), col=\"blue\", alpha=1/2)+ylim(ymin, ymax) +xlim(xmin-3, xmax)" >>script.R
echo "pdf(\"output.pdf\",width=10, height=8)">>script.R
echo "print(p+ geom_dl(aes(colour = data\$\"Category\"), method=list(\"first.bumpup\")) + scale_colour_hue(name=\"Category\") + labs(size=\"Phenotype SNPs\", color=\"Category\") + scale_size(range = c(0,50)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face=\"bold\")) + scale_color_manual(values = c(wes_palette(\"Cavalcanti1\"), wes_palette(\"Royal1\"), wes_palette(\"GrandBudapest2\"), wes_palette(\"Royal2\"), wes_palette(\"Darjeeling2\"), wes_palette(\"Zissou1\")))) ">> script.R
echo "dev.off()">>script.R


chmod 775 script.R
./script.R
rm script.R
rm *.gwascatalog* GwasCatalog.bed GwasCatalog.bed1
