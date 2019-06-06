# check read length of fastq
zcat sample.fastq.gz | awk '{if(NR%4==2) print length($1)}'


# indexing bowtie2 | bwa
bowtie2-build --threads 8 references/hg38.fa indexes/hg38 
bwa index references/hg38.fa -p indexes/hg38 


# mapping bowtie2 | bwa
bowtie2 --threads 8 -x indexes/hg38 -1 r1.fastq.gz -2 r2.fastq.gz -S sample.sam
bwa mem -t 8 indexes/hg38 r1.fastq.gz r2.fastq.gz > sample.sam


# view | sort | bam
samtools view -@ 8 -b sample.sam | samtools sort -@ 8 -O bam -o sample.sorted.bam
samtools index sample.sorted.bam


# check alignment stats
samtools view -c sample.sorted.bam
samtools flagstat sample.sorted.bam
samtools flagstat: counts for 13 categories based primarily on bit flags in the FLAG field.		
view -f	view -F	5639459 + 0 in total (QC-passed reads + QC-failed reads)
0x100	NONE	1743 + 0 secondary
0x800	NONE	0 + 0 supplementary
0x400	NONE	0 + 0 duplicates
NONE	0x4	5633869 + 0 mapped (99.90% : N/A)
0x1	0x100,400,800	5637716 + 0 paired in sequencing
0x41	0x100,400,800	2818858 + 0 read1
0x81	0x100,400,800	2818858 + 0 read2
0x2	0x100,400,800	4851102 + 0 properly paired (86.05% : N/A)
NONE	NONE	5630088 + 0 with itself and mate mapped
NONE	NONE	2038 + 0 singletons (0.04% : N/A)
NONE	NONE	71560 + 0 with mate mapped to a different chr
NONE	NONE	28422 + 0 with mate mapped to a different chr (mapQ>=5)


# bcftools mpileup | call
bcftools mpileup --threads 8 -f references/hg38.fa sample.sorted.bam -Ou | bcftools call --threads 8 -vm -Oz > sample.vcf.gz
bcftools mpileup --threads 8 -f references/hg38.fa sample.sorted.bam -Ob > sample.bcf
bcftools call --threads 8 -vm sample.bcf -Oz > sample.vcf.gz


# Remove 'chr' from chromosome names in the vcf file
zcat sample.vcf.gz | awk '{gsub(/^chr/,""); print}' | bgzip > sample.no_chr.vcf.gz


# bcftools filter
bcftools filter -i '%QUAL>=10 & INFO/DP>=10 & INFO/DP<10' sample.vcf.gz -Oz -o sample.qual_10.vcf.gz
bcftools filter -s QUAL_10 -i '%QUAL>=10' sample.vcf.gz -Oz -o sample.qual_10.vcf.gz


# bcftools summary stats
tabix -p vcf sample.vcf.gz
bcftools stats -F references/hg38.fa sample.vcf.gz > sample.vcf.gz.stats


# Compare and summarize two vcf files
bcftools stats sample_1.vcf.gz sample_2.vcf.gz > compare_vcf.sample_1_2.stats
bcftools isec -p compare.sample_1_2 sample_1.vcf.gz sample_2.vcf.gz


# PICARD AddOrReplaceReadGroups
command = "java -jar $PICARD AddOrReplaceReadGroups I={0} O={1}".format(bam, out_bam)
command += " RGID=NIPT_100"
command += " RGSM=NIPT_100"
command += " RGLB=unknown"
command += " RGPL=unknown"
command += " RGPU=unknown"
command += " VALIDATION_STRINGENCY=LENIENT"


# GATK HaplotypeCaller
gatk --java-options "-Xmx2G" HaplotypeCaller --reference /data/nh2tran/GeneSolutions/references/hg38.fa --input sample.bam --output sample.vcf.gz


