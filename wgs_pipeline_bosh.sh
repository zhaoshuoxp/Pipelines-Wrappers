picard=/software/picard-tools/1.92
gatk=/usr/bin/GenomeAnalysisTK.jar
hg19=/srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk_bundle=/srv/persistent/bliu2/shared/gatk_bundle_2.8_hg19
sample=placeholder1
reads1=placeholder2
reads2=placeholder3

# copy from /mnt/ to /srv/persistent/:
if [[ ! -d /srv/persistent/bliu2/HCASMC_eQTL/data/$sample ]]; then
	mkdir /srv/persistent/bliu2/HCASMC_eQTL/data/$sample
fi 
cp /mnt/data/WGS_HCASMC/$sample/$reads1 /srv/persistent/bliu2/HCASMC_eQTL/data/$sample/$reads1
cp /mnt/data/WGS_HCASMC/$sample/$reads2 /srv/persistent/bliu2/HCASMC_eQTL/data/$sample/$reads2

# change wd: 
cd /srv/persistent/bliu2/HCASMC_eQTL/data/$sample

# trim adapters:
cutadapt -m 30 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${reads1/gz/trimmed.gz} -p ${reads2/gz/trimmed.gz} $reads1 $reads2 
touch ${reads1/gz/trimmed.gz.done}
touch ${reads2/gz/trimmed.gz.done}

# align, sort and mark duplicates: 
bwa mem -M -t 8 -R '@RG\tID:placeholder1\tSM:placeholder1\tPL:illumina\tLB:lib1\tPU:unit1' /srv/persistent/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa ${reads1/gz/trimmed.gz} ${reads2/gz/trimmed.gz} > aln.sam 
touch aln.sam.done
java -Xmx2g -jar $picard/SortSam.jar INPUT=aln.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate 
touch sorted_reads.bam.done 
java -Xmx2g -jar $picard/MarkDuplicates.jar INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt
touch dedup_reads.bam.done 
java -Xmx2g -jar $picard/BuildBamIndex.jar INPUT=dedup_reads.bam
touch dedup_reads.bai.done 

# indel realignment: 
java -Xmx16g -jar $gatk -T RealignerTargetCreator -nt 8 -R $hg19 -I dedup_reads.bam -known $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -known $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o realignment_targets.list
touch realignment_targets.list.done
java -Xmx4g -jar $gatk -T IndelRealigner -R $hg19 -I dedup_reads.bam -targetIntervals realignment_targets.list -known $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -known $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o realigned_reads.bam 
touch realigned_reads.bam.done

# Base recalibration: 
java -Xmx4g -jar $gatk -T BaseRecalibrator -nct 8 -R $hg19 -I realigned_reads.bam -knownSites $gatk_bundle/dbsnp_138.hg19.vcf -knownSites $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -knownSites $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o recal_data.table
touch recal_data.table.done
# java -Xmx2g -jar $gatk -T BaseRecalibrator -R $hg19 -I realigned_reads.bam -knownSites $gatk_bundle/dbsnp_138.hg19.vcf -knownSites $gatk_bundle/1000G_phase1.indels.hg19.sites.vcf -knownSites $gatk_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -BQSR recal_data.table -o post_recal_data.table 
# java -Xmx2g -jar $gatk -T AnalyzeCovariates  -R $hg19 -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf
java -Xmx4g -jar $gatk -T PrintReads -nct 8 -R $hg19 -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam 
touch recal_reads.bam.done 

# Variant discovery (single sample):
## java -Xmx2g -jar $gatk -T HaplotypeCaller -R $hg19 -I recal_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o raw_variants.vcf 

# Joint variant discovery: 
java -Xmx16g -jar $gatk -T HaplotypeCaller -nct 24 -R $hg19 -I recal_reads.bam --genotyping_mode DISCOVERY --emitRefConfidence GVCF -o raw_variants.g.vcf 
touch raw_variants.g.vcf.done
