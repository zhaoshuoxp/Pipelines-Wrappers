#!/bin/bash
#####################################
# Usage:  plar.sh OUTPUTDIR labelC1,labelC2.. C1_rep1.gtf C1_rep2.gtf.. C2_rep1.gtf C2_rep2.gtf.. C1_rep1.bam,C1_rep2.bam.. C2_rep1.bam,C2_rep2.bam..      #
# Manual:                           #
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
cpc2=/home/quanyi/app/CPC2-beta/bin/CPC2.py

# generate GTF files list and BAM file args
for i in $*
do
	if [ "${i##*.}"x = "gtf"x ];then
		echo $i >> $output_dir/gtf.list
	fi
	if [ "${i##*.}"x = "bam"x ];then
		bam=${bam}" "${i}
	fi
done

# stringtie merge multi gtf files
stringtie --merge -G $plar_path/hg19_ensembl.gtf -o $output_dir/merged/merged.gtf $output_dir/gtf.list
rm $output_dir/gtf.list

# cuffdiff for fpkm 
cuffdiff -p 8 -q -L $label -o $output_dir/merged $output_dir/merged/merged.gtf $bam

#Process the Cuffmerge file: The first step traverses the cuffmerge output (merged.gtf): 
$plar plar_process_cuffmerge $output_dir/ hg19 $output_dir/merged/merged.gtf $plar_path/hg19_ensembl.genes $plar_path/hg19.size

#Annotate the transcriptome and identify potential lincRNAs:
$plar plar_identify_lincs $output_dir/ hg19 $label $plar_path/hg19_ensembl.genes $plar_path/hg19_ensembl.info1 $plar_path/hg19_ensembl.info2 $plar_path/hg19.refseq_genes $plar_path/hg19.refseq_other NONE $output_dir/merged/isoforms.fpkm_tracking $plar_path/hg19.repeat 200 2000 0.1 5 TRUE 

#Analyze the ORF sizes in candidate lincRNAs and prepare input files for CPC, HMMer
$plar plar_analyze_sequences $output_dir/ hg19 $plar_path/hg19.2bit $plar_path/hg19_mask.2bit $(pwd)/$output_dir/sequences/ $cpc2/

#Run CPĆ2
cd $output_dir/sequences
cat *.fa > all.fa
$cpc2 -i all.fa -o cpc2.result
awk '{if ($8=="coding"){print $1}}' cpc2.result |cut -f 3 -d '|' > ../coding.list

#Run HMMER
hmmsearch --cpu 32 --noali --tblout hmm.result -E 0.001 -o hmm.log $plar_path/Pfam/Pfam-A.hmm all.fa
cut -f 3 -d '|' hmm.result |grep -v ^# >> ../coding.list
cd ..
rm -r sequences

#Remove protein coding potential lncRNA (combined HMMER and CPC2 results)
gunzip -c hg19.lincs.f1.bed.gz |grep -v -f coding.list |gzip > hg19.lincs.f1.clean.bed.gz
rm coding.list
cd ..

#Further filter transcripts near coding genes etc.
$plar plar_filter $output_dir/ hg19 $plar_path/hg19_ensembl.genes $plar_path/hg19_ensembl.info1 $plar_path/hg19_ensembl.info2 $plar_path/hg19.gap $plar_path/hg19.size 500 2000 TRUE

rm cpc_names.txt

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################