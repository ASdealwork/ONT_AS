#!/bin/bash

#下机数据做basecall
#Perform basecalling on raw sequencing data
dorado basecaller /home/xiary/my_data/software/dorado-0.5.1-linux-x64/model/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
	/data/ONT_data/example_sample/20231207_1555_6E_PAQ53476_984654b1/pod5 --emit-fastq --recursive -x 'cuda:all' \
	> /data/ONT_data/example_sample/20231207_1555_6E_PAQ53476_984654b1/basecall/example_sample.fastq

#NanoFilt过滤Q-score <10的reads
#Filter reads with Q-score < 10 using NanoFilt
NanoFilt -q 10 /data/ONT_data/example_sample/20231207_1555_6E_PAQ53476_984654b1/basecall/example_sample.fastq >/data/ONT_data/example_sample/analysis/example_sample_q10.fastq

#minimap2做比对
#Align reads using minimap2
minimap2 -t 50 -I 16g -ax map-ont /data/ONT_data/example_sample/ref/file1.Col-PEK1.5_Chr1-5_20220523.fasta /data/ONT_data/example_sample/analysis/example_sample_q10.fastq > /data/ONT_data/example_sample/analysis/example_sample_q10.sam

grep '^@' /data/ONT_data/example_sample/analysis/example_sample_q10.sam > /data/ONT_data/example_sample/analysis/example_sample_q10_header.sam

awk '{if (($2 == 0 || $2 == 16) && $5 > 30) print}' /data/ONT_data/example_sample/analysis/example_sample_q10_header.sam > /data/ONT_data/example_sample/analysis/example_sample_q10_header_body.sam

cat /data/ONT_data/example_sample/analysis/example_sample_q10_header.sam /data/ONT_data/example_sample/analysis/example_sample_q10_header_body.sam > /data/ONT_data/example_sample/analysis/example_sample_q10_filter.sam

rm -r /data/ONT_data/example_sample/analysis/example_sample_q10_header.sam /data/ONT_data/example_sample/analysis/example_sample_q10_header_body.sam

samtools view -@ 50 -hSb example_sample_q10_filter.sam | samtools sort -@ 50 - > example_sample_q10_filter_sorted.bam

samtools index -c example_sample_q10_filter_sorted.bam

#STAM-seq输入文件预处理
#Preprocess STAM-seq input files
##取出bam文件中某段区域的reads比对情况
#Extract alignment information of reads from a specific region in the BAM file
samtools view -b -L regions.bed example_sample_q10_filter_sorted.bam > example_sample_q10_filter_sorted_target.bam
samtools view -b -U example_sample_q10_filter_sorted_nontarget.bam -L regions.bed example_sample_q10_filter_sorted.bam
samtools index -c example_sample_q10_filter_sorted_target.bam
samtools index -c example_sample_q10_filter_sorted_nontarget.bam

samtools view example_sample_q10_filter_sorted_target.bam| awk '{print length($10)}' > target_length.txt
samtools view example_sample_q10_filter_sorted_nontarget.bam | awk '{print length($10)}' > nontarget_length.txt

##统计reads覆盖数目及长度
##Count the number and length of reads coverage
bamCoverage --bam example_sample_q10_filter_sorted.bam --outFileName example_sample_q10_filter_sorted.bw --binSize 1 --normalizeUsing None --skipNonCoveredRegions --numberOfProcessors 16

#端粒区域
#Telomere region
##获取比对到指定区域的reads
##Retrieve reads aligned to the specified region
samtools view -b example_sample_q10_filter_sorted.bam Chr04:1-150000 > example_sample_q10_filter_sorted_target.bam
##获取指定ID的序列
##Retrieve sequences with specified IDs
samtools view example_sample_q10_filter_sorted_target.bam |cut -f1|sort|uniq example_sample_align_T.id.txt
cat example_sample_align_T.id.txt |while read id;do grep -A 1 ${id} /data/ONT_data/example_sample/analysis/example_sample_q10.fastq >> example_sample_align_T.reads.fa;done

sed -i 's#^@#>#g' example_sample_align_T.reads.fa
##统计比对到端粒区域reads含有的CCCTAAA数目
##Count the number of reads aligned to the telomere region containing the sequence CCCTAAA
python CCCTAAA_stat.num.py example_sample_align_T.reads.fa >CCCTAAA_stat.num_out.txt


#统计fastq.gz文件中reads长度信息,将结果文件作为cal_length_data.py的输入，进而计算测序reads长度分布信息
#Calculate the length distribution of sequencing reads by extracting length information from fastq.gz files, and then use the resulting file as input for cal_length_data.py
seqkit fx2tab -j 30 -l  -n -i -H /data/ONT_data/example_sample/20231207_1555_6E_PAQ53476_984654b1/basecall/example_sample.fastq > example_sample.reads.sequence_length_num.txt

