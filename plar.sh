#!/bin/bash
#####################################
# Usage:  plar.sh OUTPUTDIR labelC1,labelC2.. fr|rf|un(null) C1_rep1.gtf C1_rep2.gtf.. C2_rep1.gtf C2_rep2.gtf.. C1_rep1.bam,C1_rep2.bam.. C2_rep1.bam,C2_rep2.bam..      #
#####################################
# setup PATH
output_dir=$1
label=$2

if [ ! -d $output_dir ];then 
mkdir -p $1/sequences
fi
current=$(pwd)
plar_path=/home/quanyi/app/PLAR
plar=$plar_path/find_lincs.csh
cpc2_path=/home/quanyi/app/CPC2-beta
cpc2=$cpc2_path/bin/CPC2.py

# generate GTF files list and BAM file args
for i in $*
do
	if [ "${i##*.}"x = "gtf"x ];then
		echo $i >> $output_dir/gtf.list
	elif [ "${i##*.}"x = "bam"x ];then
		bam=${bam}" "${i}
	fi
done

# deternmine the strand direction
if [ $3 = 'fr' ];then
	strand='fr-secondstrand'
elif [ $3 = 'rf' ];then
	strand='fr-firststrand'
elif [ $3 = 'un' ];then
	strand='fr-unstranded'
else
	strand='fr-unstranded'
fi

# stringtie merge multi gtf files
stringtie --merge -G $plar_path/hg19_ensembl.gtf -F 0 -T 0 -f 0 -c 0 -o $output_dir/merged/merged.gtf $output_dir/gtf.list
rm $output_dir/gtf.list

# cuffdiff for fpkm 
cuffdiff -p 8 -q -L $label -o $output_dir/merged --library-typ $strand $output_dir/merged/merged.gtf $bam

#Process the Cuffmerge file: The first step traverses the cuffmerge output (merged.gtf): 
$plar plar_process_cuffmerge $output_dir/ hg19 $output_dir/merged/merged.gtf $plar_path/hg19_ensembl.genes $plar_path/hg19.size

#Annotate the transcriptome and identify potential lincRNAs:
$plar plar_identify_lincs $output_dir/ hg19 $label $plar_path/hg19_ensembl.genes $plar_path/hg19_ensembl.info1 $plar_path/hg19_ensembl.info2 $plar_path/hg19.refseq_genes $plar_path/hg19.refseq_other NONE $output_dir/merged/isoforms.fpkm_tracking $plar_path/hg19.repeat 200 200 0.1 1 TRUE 

#Analyze the ORF sizes in candidate lincRNAs and prepare input files for CPC, HMMer
$plar plar_analyze_sequences $output_dir/ hg19 $plar_path/hg19.2bit $plar_path/hg19_mask.2bit $(pwd)/$output_dir/sequences/ $cpc2_path/

# custom code for CPC2 and HMMER
##########################
#Run CPC2
#cd $output_dir/sequences
#$cpc2 -i all.fa -o cpc2.result
#awk '{if ($8=="coding"){print $1}}' cpc2.result |cut -f 3 -d '|' > ../coding.list

#Run HMMER
#cat *.faa >all.faa
#hmmscan --cpu 16 --noali --tblout hmm.result -E 0.001 -o hmm.log $plar_path/Pfam/Pfam-A.hmm all.faa
#cut -f 3 -d '|' hmm.result |grep -v ^# >> ../coding.list

#Remove protein coding potential lncRNA (combined HMMER and CPC2 results)
#gunzip -c hg19.lincs.f1.bed.gz |grep -v -f coding.list |gzip > hg19.lincs.f1.clean.bed.gz
#rm coding.list
#cd ../..
#########################

cd $output_dir/sequences
cat *.fa > all.fa
$cpc2 -i all.fa -o cpc2.result
awk -v OFS="\t" '{print $1,$2,$8,$4}' cpc2.result > cpc.result
rm $cpc2_path/run_hg19.csh

#Run HMMER
cat *.faa >all.faa
hmmscan --cpu 16 -o hmmer.log --noali -E 0.001 --tblout hmmer.result $plar_path/Pfam/Pfam-A.hmm all.faa
rm hg19_hmmer_all.csh
cd ../..

#Combine the results of the coding predictors
$plar plar_parse_predictors $output_dir/ hg19 NONE $output_dir/sequences/hmmer.result $output_dir/sequences/cpc.result

#Further filter transcripts near coding genes etc.
$plar plar_filter $output_dir/ hg19 $plar_path/hg19_ensembl.genes $plar_path/hg19_ensembl.info1 $plar_path/hg19_ensembl.info2 $plar_path/hg19.gap $plar_path/hg19.size 500 2000 TRUE

#Add single-exon lncRNA back
gunzip -c $output_dir/hg19.lincs.f1.pure.bed.gz > $output_dir/hg19.lincs.f1.pure.bed
gunzip -c $output_dir/hg19.lincs.f1.with_filters.bed.gz > $output_dir/hg19.lincs.f1.with_filters.bed
awk '$10==1' $output_dir/hg19.lincs.f1.with_filters.bed >> $output_dir/hg19.lincs.f1.pure.bed
sort -u $output_dir/hg19.lincs.f1.pure.bed | sort -k1,1 -k2,2n > final_lncRNA.bed

rm cpc_names.txt

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################