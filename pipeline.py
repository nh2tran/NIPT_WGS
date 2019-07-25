from ngs_modules import *

#mapping_bwa(index, fastq_dir, output_dir, read_group_dict, num_threads=1, num_test=1)


# to include into mapping_bwa???
#process_multi_samples(gatk_mark_duplicates, bam_list_file, output_dir, num_threads=1, num_test=1)


# to include into mapping_bwa???
#process_multi_samples(gatk_build_bam_index, bam_list_file, output_dir, num_threads=1, num_test=1)


# to include into mapping_bwa???
# samtools view -q 10/30???


#samtools_flagstat(bam_list_file, output_dir, num_threads=1, num_test=1)


# to include into mapping_bwa???
# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp/filter_2plus_reads/"
# ~ num_threads = 16
# ~ num_test = None
# ~ process_multi_samples(filter_2plus_reads, bam_list_file, output_dir, num_threads, num_test)


# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp/filter_2plus_reads.list"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp/filter_2plus_reads_hist_genome/"
# ~ num_threads = 16
# ~ num_test = None
# ~ genomecov_hist(bam_list_file, output_dir, num_threads, num_test)


# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
# ~ ref_dict_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.dict"
# ~ ref_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
# ~ dbSNP_vcf = "/data/nh2tran/GeneSolutions/dbSNP/human_9606_b150_GRCh38p7/All_20170710.with_chr.vcf.gz"
# ~ output_dir = "/data/nh2tran/GeneSolutions/temp/gatk_bqsr_workflow/"
# ~ num_threads = 16
# ~ gatk_bqsr_workflow(bam_list_file, ref_dict_file, ref_file, dbSNP_vcf, output_dir, num_threads)


# ~ bam_list_file = "/data/nh2tran/GeneSolutions/temp/filter_2plus_reads.list"
# ~ reference_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
# ~ output_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.filter_2plus_reads.Mutect2.vcf.gz"
# ~ gatk_multi_chromosome(bam_list_file, reference_file, output_file, num_threads=8)


# ~ input_1000genomes = "/data/nh2tran/GeneSolutions/1000genomes/release_20130502_GRCh38_positions/"
# ~ population = "khv"
# ~ sample_file = "/data/nh2tran/GeneSolutions/1000genomes/khv/khv_samples.txt"
# ~ output_dir = "/data/nh2tran/GeneSolutions/1000genomes/khv/"
# ~ num_threads = 24
# ~ extract_sample_variants(input_1000genomes, population, sample_file, output_dir, num_threads)


print("Pipeline executed successfully.")

