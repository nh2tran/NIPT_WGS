from ngs_modules import *

########################################################################
# step 1 - preprocessing: fastqc | trim | align | sort | dedup | index
# ~ index = "/data/nh2tran/GeneSolutions/indexes/hg38_selected"
# ~ fastq_dir = "/data/nh2tran/GeneSolutions/fastq/NIPT_683/"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_683.step_1_preprocessing/"
# ~ read_group_dict = {'SM': "NIPT_700", 'PL': "unknown", 'PU': "unknown"}
# ~ num_samples = 6 # 4 threads, 16G
# ~ num_test = None
# ~ #preprocessing(index, fastq_dir, output_dir, read_group_dict, num_samples, num_test)
########################################################################
# summarize flagstat of mapped reads for raw
# ~ num_samples = 6 # 4 threads
# ~ num_test = None
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_1_preprocessing.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_5_flagstat_raw/"
# ~ #samtools_flagstat(bam_list_file, output_dir, num_samples, num_test)


########################################################################
# step 2 - genomecov_hist: summarize percentage of nonzero, 1x, 2x, etc, and depth per bam.
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_1_preprocessing.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_2_genomecov_hist/"
# ~ num_samples = 24 # 1 threads
# ~ num_test = None
# ~ #genomecov_hist(bam_list_file, output_dir, num_samples, num_test)
########################################################################


########################################################################
# step 2 - genomecov_2x:  summarize distribution of 2x regions across genome and samples.
# ~ hg38_selected_genome = "/data/nh2tran/GeneSolutions/references/hg38_selected.genome"
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_1_preprocessing.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_2_genomecov_2x/"
# ~ num_samples = 24 # 1 threads
# ~ num_test = None
# ~ #genomecov_2x(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test)
########################################################################


########################################################################
# ~ # step 3 - filter_q30_1read:  filter alignments with mapq >= 30 and sample 1 overlapping read.
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_1_preprocessing.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_3_filter_q30_1read/"
# ~ num_samples = 24 # 1 threads
# ~ num_test = None
# ~ os.system("mkdir " + output_dir)
# ~ #process_multi_samples(filter_q30_1read, bam_list_file, output_dir, num_samples, num_test)
########################################################################


########################################################################
# step 4 - genomecov_bg: calculate aggregated bed graph and sequencing depth
# ~ hg38_selected_genome = "/data/nh2tran/GeneSolutions/references/hg38_selected.genome"
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_3_filter_q30_1read.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_4_genomecov_bg/"
# ~ num_samples = 24 # 1 threads
# ~ num_test = None
# ~ #genomecov_bg(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test)
########################################################################


########################################################################
# step 5 - summary
########################################################################
# summarize and draw distribution of sequencing depth
# ~ hist_file = "/data/nh2tran/GeneSolutions/temp/NIPT_2683.step_4_genomecov_bg/merged.genomecov_hist.full"
# ~ #genomecov_hist_sum(hist_file)
########################################################################
# summarize flagstat of mapped reads for filter
# ~ num_samples = 6 # 4 threads
# ~ num_test = None
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_3_filter_q30_1read.bam_list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_5_flagstat_filter/"
# ~ #samtools_flagstat(bam_list_file, output_dir, num_samples, num_test)
########################################################################
# igv tracks: depth, mapq (NIPT_100.q_10), k50.umap
# "ls *.bw > temp.merged_chr_bw.list"
# "bigWigMerge -inList temp.merged_chr_bw.list temp.merged.bg"
# "mkdir temp"
# "sort -k1,1 -k2,2n -T temp --parallel 16 temp.merged.bg > temp.merged.sorted.bg"
# "bedGraphToBigWig temp.merged.sorted.bg /data/nh2tran/GeneSolutions/references/hg38_selected.genome merged.bw"


########################################################################
# step n - Mutect2: variant calling
########################################################################
# Note that gatk requires the bam_list file to have extension ".list"
# ~ reference_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_3_filter_q30_1read.bam_list.list"
# ~ output_file = "/data/nh2tran/GeneSolutions/temp_1K/NIPT_2683.step_n_Mutect2.vcf.gz"
# ~ num_threads = 10 # 8G
# need to merge Mutect2 stats files, FilterMutectCalls, and annotate NIPT_AF
# ~ #gatk_multi_chromosome(bam_list_file, reference_file, output_file, num_threads)
########################################################################
# step n - FilterMutectCalls
# ~ gatk FilterMutectCalls --variant test.vcf.gz --reference ../references/hg38_selected.fa --output test.filtered.vcf.gz --stats test.vcf.gz.merged.stats
########################################################################
# step n - annotate NIPT_AF
# ~ bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' NIPT_2683.step_n_FilterMutectCalls.vcf.gz | bgzip > temp_annotate.tab.gz
# ~ tabix -s1 -b2 -e2 temp_annotate.tab.gz
# ~ vim temp_annotate.hdr
# ~ ##INFO=<ID=NIPT_AF,Number=A,Type=Float,Description="Estimated allele frequency from NIPT, range (0,1)">
# ~ bcftools annotate -a temp_annotate.tab.gz -h temp_annotate.hdr -c CHROM,POS,REF,ALT,NIPT_AF NIPT_2683.step_n_FilterMutectCalls.vcf.gz -Oz -o NIPT_2683.step_n_NIPT_AF.vcf.gz
# ~ rm temp_annotate.*
########################################################################


########################################################################
# step z - filter variants
#bcftools view --include 'INFO/DP>=162 & INFO/DP<=566' NIPT_2683.step_n_NIPT_AF.vcf.gz
# ~ bcftools view NIPT_2683.step_n_NIPT_AF.vcf.gz --targets ^chrY --exclude 'FILTER="weak_evidence" || FILTER="strand_bias" || FILTER="contamination"' -Oz -o NIPT_2683.vcf.gz
# ~ tabix -p vcf NIPT_2683.vcf.gz
# ~ bcftools stats NIPT_2683.vcf.gz > NIPT_2683.vcf.gz.stats
# step z - filter SNPs, LeftAlignAndTrimVariants
# ~ bcftools norm --multiallelics - --fasta-ref ../references/hg38_selected.fa NIPT_2683.vcf.gz | bcftools view --types snps -Oz -o NIPT_2683.snps_only.vcf.gz
# ~ tabix -p vcf NIPT_2683.snps_only.vcf.gz
# ~ bcftools stats NIPT_2683.snps_only.vcf.gz > NIPT_2683.snps_only.vcf.gz.stats
########################################################################


########################################################################
# Analysis of variant calling (Figure 2)
# extract KHV, EAS variants from 1000genomes
# ~ input_1000genomes = "/data/nh2tran/GeneSolutions/1000genomes/release_20130502_GRCh38_positions/"
# ~ population = "all"
# ~ sample_file = "/data/nh2tran/GeneSolutions/1000genomes/khv/khv_samples.txt"
# ~ output_dir = "/data/nh2tran/GeneSolutions/1000genomes/khv/"
# ~ num_threads = 24
# ~ #extract_sample_variants(input_1000genomes, population, sample_file, output_dir, num_threads)


# add 'chr' and add contig length
# ~ zcat khv.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | bgzip > khv.temp_1.vcf.gz
# ~ gatk UpdateVcfSequenceDictionary --INPUT khv.temp_1.vcf.gz --SEQUENCE_DICTIONARY ../references/hg38_selected.dict --OUTPUT khv.vcf
# ~ bgzip khv.vcf
# ~ tabix -p vcf khv.vcf.gz
# ~ bcftools stats khv.vcf.gz > khv.vcf.gz.stats
# ~ rm khv.temp_1.vcf.gz
# filter SNPs, LeftAlignAndTrimVariants
# ~ bcftools norm --multiallelics - --fasta-ref ../references/hg38_selected.fa khv.vcf.gz | bcftools view --types snps -Oz -o khv.snps_only.vcf.gz
# ~ tabix -p vcf khv.snps_only.vcf.gz
# ~ bcftools stats khv.snps_only.vcf.gz > khv.snps_only.vcf.gz.stats

# annotate KHV, EAS, dbSNPS
# ~ gatk --java-options '-Xmx32G' VariantAnnotator --variant NIPT_2683.snps_only.vcf.gz --dbsnp dbSNP_151.vcf.gz --output NIPT_2683.snps_only.annotated.temp_1.vcf.gz
# ~ gatk --java-options '-Xmx32G' VariantAnnotator --variant NIPT_2683.snps_only.dbSNP_151.vcf.gz --comp:1KG_KHV khv.snps_only.vcf.gz --comp:1KG_EAS eas.snps_only.vcf.gz --output NIPT_2683.snps_only.annotated.vcf.gz
# ~ rm NIPT_2683.snps_only.annotated.temp_1.vcf.gz*


# compare NIPT_2683 to KHV, EAS, dbSNPS
# ~ bcftools stats NIPT_2683.snps_only.vcf.gz khv.snps_only.vcf.gz > compare_snps.NIPT_2683.khv.stats
# ~ draw_figure_2_venn()


# AF distribution with respect to Venn diagram
# ~ bcftools view NIPT_2683.snps_only.annotated.vcf.gz --include 'INFO/1KG_KHV == 1' -Oz -o NIPT_2683.snps_only.venn_khv1.vcf.gz
# ~ bcftools view NIPT_2683.snps_only.annotated.vcf.gz --include 'INFO/1KG_KHV == 0 && (INFO/1KG_EAS == 1 || ID != ".")' -Oz -o NIPT_2683.snps_only.venn_khv0_eas1_or_dbSNP1.vcf.gz
# ~ bcftools view NIPT_2683.snps_only.annotated.vcf.gz --include 'INFO/1KG_KHV == 0 && INFO/1KG_EAS == 0 && ID == "."' -Oz -o NIPT_2683.snps_only.venn_novel.vcf.gz
# ~ bcftools stats NIPT_2683.snps_only.venn_novel.vcf.gz --af-tag NIPT_AF --af-bins af_bins.txt > NIPT_2683.snps_only.venn_novel.vcf.gz.stats


# compare allele frequency between NIPT and KHV 
# ~ bcftools isec NIPT_2683.snps_only.vcf.gz khv.snps_only.vcf.gz -p isec_snps.NIPT_2683.khv
# ~ isec_dir = "/data/nh2tran/GeneSolutions/temp/isec_snps.NIPT_2683.khv/"
# ~ #sample_AF_file, population_AF_file = get_AF(isec_dir)
# ~ sample_AF_file = isec_dir + "sample_AF_list.npy"
# ~ population_AF_file = isec_dir + "population_AF_list.npy"
# ~ AF_max = 1.0
# ~ scatter_plot_res = 100
# ~ compare_AF(sample_AF_file, population_AF_file, AF_max, scatter_plot_res)


########################################################################


print("Pipeline executed successfully.")
sys.exit()


# ~ bed_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.bedgraph"
# ~ hist_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.hist"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp/create_chr_mask/"
# ~ #create_chr_mask(bed_file, hist_file, output_dir)


# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
# ~ ref_dict_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.dict"
# ~ ref_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
# ~ dbSNP_vcf = "/data/nh2tran/GeneSolutions/dbSNP/human_9606_b150_GRCh38p7/All_20170710.with_chr.vcf.gz"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp/gatk_bqsr_workflow/"
# ~ num_threads = 16
# ~ #gatk_bqsr_workflow(bam_list_file, ref_dict_file, ref_file, dbSNP_vcf, output_dir, num_threads)


