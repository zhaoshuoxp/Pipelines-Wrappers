#!/bin/bash

#create bed from bam, requires bedtools bamToBed
bamToBed -i $1 -split > accepted_hits.bed

#create plus and minus strand bedgraph
cat accepted_hits.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=$len_hg19 stdin | sort -k1,1 -k2,2n > accepted_hits.bedGraph

bedGraphToBigWig accepted_hits.bedGraph $len_hg19 $1.bw

#removing intermediery files
rm accepted_hits.bed
rm accepted_hits.bedGraph