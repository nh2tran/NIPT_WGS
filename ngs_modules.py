from __future__ import print_function


import os
import re
import math
import numpy as np
import pandas as pd
import csv
import multiprocessing

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
matplotlib.rcParams.update({'font.size': 11})

from Bio import SeqIO
from Bio.SeqIO import FastaIO


NUM_THREADS = 16
hg38_fasta = "/data/nh2tran/GeneSolutions/references/hg38.fa"


def compute_allele_frequency():
  #os.system("bcftools view isec/0002.vcf -H | cut -f 1,2,10 > isec/gatk_AD.txt")
  #os.system("bcftools view isec/0003.vcf -H | cut -f 8 > isec/khv_AF.txt")
  #os.system("paste -d '\t' isec/gatk_AD.txt isec/khv_AF.txt > isec/gatk_AD_khv_AF.txt")
  # ~ AD_freq_list = []
  # ~ AF_freq_list = []
  # ~ skip = 0
  # ~ with open("isec/gatk_AD_khv_AF.txt", 'r') as handle:
    # ~ for line in handle.readlines(): # chr | position | gatk | khv
      # ~ AD = line.split('\t')[2].split(':')[1] # ref_count, alt_count, alt_count
      # ~ AF = [x.split('=')[1] for x in line.split('\t')[3].split(';') if x[:3] == 'AF='] # alt_freq, alt_freq
      # ~ assert len(AF) == 1, "Error: wrong AF"
      # ~ AF = AF[0]
      # ~ AD_count = AD.split(',')
      # ~ AF_freq = AF.split(',')
      # ~ assert len(AD_count) - 1 == len(AF_freq), "Error: AD and AF not matched"
      # ~ AD_count = [int(x) for x in AD_count]
      # ~ total_count = float(sum(AD_count))
      # ~ if total_count == 0:
        # ~ skip += 1
        # ~ continue
      # ~ AF_freq = [float(x) for x in AF_freq]
      # ~ for alt_count, alt_freq in zip(AD_count[1:], AF_freq):
        # ~ AD_freq_list.append(alt_count / total_count)
        # ~ AF_freq_list.append(alt_freq)
  # ~ print("skip =", skip)
  # ~ print("len(AD_freq_list) =", len(AD_freq_list))
  # ~ print("len(AF_freq_list) =", len(AF_freq_list))
  # ~ np.save("AD_freq_list.npy", AD_freq_list)
  # ~ np.save("AF_freq_list.npy", AF_freq_list)
  AD_freq_list = np.load("AD_freq_list.npy")
  AF_freq_list = np.load("AF_freq_list.npy")
  print("AF_freq_list.shape =", AF_freq_list.shape)
  corr = np.corrcoef(AD_freq_list, AF_freq_list)[0, 1]
  mse = np.median((AD_freq_list - AF_freq_list) ** 2)
  mae = np.median(np.absolute(AD_freq_list - AF_freq_list))
  mape = np.median(np.absolute(AD_freq_list - AF_freq_list) / ((AD_freq_list + AF_freq_list) / 2))
  print("corr =", corr)
  print("mse =", mse)
  print("mae =", mae)
  print("mape =", mape)
  down_sampling = np.random.choice(AF_freq_list.shape[0], size=AF_freq_list.shape[0]/10, replace=False)
  fig, ax = pyplot.subplots(figsize=(8, 8))
  pyplot.scatter(x=AF_freq_list[down_sampling], y=AD_freq_list[down_sampling],
                 marker='.', s=0.1)
  ax.set_xlabel('KHV AF frequency')
  ax.set_ylabel('gatk AD frequency')
  pyplot.savefig("scatter_plot.png")


# ~ total = 601102.0 + 4098436.0 + 4534279.0
# ~ venn_plot = venn2(subsets=(601102, 4098436, 4534279),
                  # ~ set_labels=['GATK', 'KHV >= 2%'],
                  # ~ subset_label_formatter=lambda x: "{0:,}\n({1:.0%})".format(x, x/total))
# ~ pyplot.savefig("venn2.png")




def gdrive_download(input_file, output_dir):
  """Download multiple files from google drive with link IDs.

     Usage:
       input_file = "googledrive_second_600_samples_754_link_IDs.txt"
       output_dir = "/data/nh2tran/GeneSolutions/fastq/temp"
       gdrive_download(input_file, output_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("gdrive_download()")
  print("input_file =", input_file)
  print("output_dir =", output_dir)
  print("".join(["="] * 80)) # section-separating line

  with open(input_file, 'r') as file_handle:
    for line in file_handle.readlines():
      link_ID = line.strip()
      command = ["gdrive download"]
      command += [link_ID]
      command += ["--path", output_dir]
      command += ["--force"]
      command = " ".join(command)
      print("command =", command)
      while True:
        status = os.system(command)
        if status == 0:
          break


def mapping_bwa(index, fastq_dir, output_dir, read_group_dict, num_threads=NUM_THREADS, num_test=None):
  """Map paired-end fastq files using bwa.

     Usage:
       index = "/data/nh2tran/GeneSolutions/indexes/hg38"
       fastq_dir = "/data/nh2tran/GeneSolutions/fastq/NIPT_600/"
       output_dir = "/data/nh2tran/GeneSolutions/temp/"
       read_group_dict = {'LB': "unknown", 'PL': "unknown", 'SM': "NIPT_600", 'PU': "unknown"}
       mapping_bwa(index, fastq_dir, output_dir, read_group_dict)
  """

  print("".join(["="] * 80)) # section-separating line
  print("Mapping with bwa")
  print("index =", index)
  print("fastq_dir =", fastq_dir)
  print("output_dir =", output_dir)
  print("read_group_dict =", read_group_dict)
  print("num_threads =", num_threads)
  print("".join(["="] * 80)) # section-separating line

  # get pairs of fastq files
  fastq_list = [os.path.join(fastq_dir, x) for x in os.listdir(fastq_dir)]
  r1_list = [x for x in fastq_list if 'R1' in x]
  r2_list = [x for x in fastq_list if 'R2' in x]
  # check paired-end matching
  r1_list.sort()
  r2_list.sort()
  if num_test is not None:
    r1_list = r1_list[:num_test]
    r2_list = r2_list[:num_test]
  for r1, r2 in zip(r1_list, r2_list):
    assert r1 == r2.replace('_R2_', '_R1_'), "Error: R1 and R2 not matched"
  print("Number of fastq pairs:", len(r1_list))
  print()

  # trim fastq; map paired reads; then sort and index bam file
  # unpaired reads are kept, paired are removed after mapping
  num_threads = str(num_threads)
  bam_dir = output_dir + "bam/"
  trim_unpaired_dir = output_dir + "trim_unpaired/"
  os.system("mkdir " + trim_unpaired_dir)
  os.system("mkdir " + bam_dir)
  for r1, r2 in zip(r1_list, r2_list):

    sample_id = r1.split('/')[-1].split('_')[0]

    trim_r1_paired = sample_id + ".R1.paired.fastq.gz"
    trim_r1_unpaired = sample_id + ".R1.unpaired.fastq.gz"
    trim_r2_paired = sample_id + ".R2.paired.fastq.gz"
    trim_r2_unpaired = sample_id + ".R2.unpaired.fastq.gz"
    # keep unpaired in trim_unpaired_dir after paired-end mapping
    trim_r1_unpaired = trim_unpaired_dir + trim_r1_unpaired
    trim_r2_unpaired = trim_unpaired_dir + trim_r2_unpaired
    # CAREFUL!!! trim_r1_paired and trim_r2_paired are removed after alignment

    command_trim = ["java -jar $trimmomatic PE"]
    command_trim += [r1, r2]
    command_trim += [trim_r1_paired, trim_r1_unpaired, trim_r2_paired, trim_r2_unpaired]
    command_trim += ["ILLUMINACLIP:/data/nh2tran/GeneSolutions/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10"]
    command_trim += ["LEADING:3"]
    command_trim += ["TRAILING:3"]
    command_trim += ["SLIDINGWINDOW:4:15"]
    command_trim += ["CROP:75"]
    command_trim += ["MINLEN:36"]
    command_trim += ["-threads", num_threads]
    command_trim = " ".join(command_trim)
    print("command_trim =", command_trim)
    print()
    # CAREFUL!!! trim_r1_paired and trim_r2_paired are removed after alignment
    os.system(command_trim)
    
    bam_file = sample_id + ".sorted.bam"
    bam_file = bam_dir + bam_file
    read_group = ["@RG"]
    read_group.append("ID:" + sample_id)
    read_group.append("LB:" + sample_id)
    read_group.append("SM:" + read_group_dict['SM'])
    read_group.append("PL:" + read_group_dict['PL'])
    read_group.append("PU:" + read_group_dict['PU'])
    read_group = "\\t".join(read_group)
    read_group = "\'" + read_group + "\'"

    command_bwa = ["bwa mem"]
    command_bwa += [index, trim_r1_paired, trim_r2_paired]
    command_bwa += ["-M"]
    command_bwa += ["-R", read_group]
    command_bwa += ["-t", num_threads]
    command_bwa = " ".join(command_bwa)

    command_view = ["samtools view"]
    command_view += ["-b"]
    command_view += ["-@", num_threads]
    command_view = " ".join(command_view)

    command_sort = ["samtools sort"]
    command_sort += [">", bam_file]
    command_sort += ["-@", num_threads]
    command_sort = " ".join(command_sort)

    command_map = " | ".join([command_bwa, command_view, command_sort])
    print("command_map =", command_map)
    print()
    os.system(command_map)

    command_index = ["samtools index"]
    command_index += [bam_file]
    command_index = " ".join(command_index)
    print("command_index =", command_index)
    print()
    os.system(command_index)

    # remove trim paired fastq after after paired-end mapping
    os.system("rm " + trim_r1_paired)
    os.system("rm " + trim_r2_paired)


def process_multi_samples(function, bam_list_file, output_dir, num_threads=1, num_test=None):
  """Process multiple samples in parallel threads.

     Usage:
       function = samtools_flagstat_per_bam
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam_list.txt"
       output_dir = "/data/nh2tran/GeneSolutions/temp/samtools_flagstat/"
       process_multi_samples(function, bam_list_file, output_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("process_multi_samples()")
  print("function =", function.__name__)
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("".join(["="] * 80)) # section-separating line

  with open(bam_list_file, 'r') as file_handle:
    bam_list = [line.strip() for line in file_handle.readlines()]
  print("Number of bam files:", len(bam_list))
  # limite the number of bam files in test mode
  if num_test is not None:
    bam_list = bam_list[:num_test]

  # create a folder to store results
  if output_dir is not None:
    assert not os.path.exists(output_dir), "output_dir exists " + output_dir
    os.makedirs(output_dir)

  output_dir_list = [output_dir] * len(bam_list)
  argument_list = zip(bam_list, output_dir_list)
  pool = multiprocessing.Pool(processes=num_threads)
  pool.map(function, argument_list)
  pool.close()
  pool.join()


def samtools_flagstat_per_bam(io_tuple):

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]
  log_file = output_dir + bam_file_name + ".flagstat"
  command = ["samtools flagstat"]
  command += [bam_file]
  command += [">", log_file]
  command = " ".join(command)
  os.system(command)


def samtools_flagstat(bam_list_file, output_dir, num_threads=1, num_test=None):
  """Summarize bam files.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam_list.txt"
       output_dir = "/data/nh2tran/GeneSolutions/temp/samtools_flagstat/"
       samtools_flagstat(bam_list_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("samtools_flagstat()")
  print("bam_list_file =", bam_list_file)
  print("".join(["="] * 80)) # section-separating line

  process_multi_samples(samtools_flagstat_per_bam, bam_list_file, output_dir, num_threads=num_threads, num_test=num_test)

  log_file_list = [output_dir + x for x in os.listdir(output_dir)]
  total_count_list = []
  for log_file in log_file_list:
    count_pair_list = []
    with open(log_file, 'r') as file_handle:
      for line in file_handle.readlines():
        line_split = line.split(' ')
        count_pair_list.append([int(line_split[0]), int(line_split[2])])
    total_count_list.append(count_pair_list)
  # sum up counts over samples
  total_count_array = np.array(total_count_list)
  total_count_array = np.sum(total_count_array, axis=0)

  print("num_total =", total_count_array[0])
  print("num_secondary =", total_count_array[1])
  print("num_supplementary =", total_count_array[2])
  print("num_duplicate =", total_count_array[3])
  print("num_mapped =", total_count_array[4])
  print("num_paired =", total_count_array[5])
  print("num_read1 =", total_count_array[6])
  print("num_read2 =", total_count_array[7])
  print("num_properly_paired =", total_count_array[8])
  print("num_itself_mate_mapped =", total_count_array[9])
  print("num_singleton =", total_count_array[10])
  print("num_mate_mapped_diff_chr =", total_count_array[11])
  print("num_mate_mapped_diff_chr_mapq_5 =", total_count_array[12])

  mapped_rate = float(total_count_array[4][0]) / total_count_array[0][0]
  properly_paired_rate = float(total_count_array[8][0]) / total_count_array[5][0]
  duplicate_rate = float(total_count_array[3][0]) / total_count_array[4][0]
  print("mapped_rate =", mapped_rate)
  print("properly_paired_rate =", properly_paired_rate)
  print("duplicate_rate =", duplicate_rate)


def gatk_mark_duplicates(io_tuple):

  bam_file, output_dir = io_tuple

  # gatk MarkDuplicates -I 10-W3790.sorted.bam -O 10-W3790.sorted.marked_duplicates.bam -M metrics.txt
  bam_file_name = bam_file.split('/')[-1]
  prefix = bam_file_name[:-4]
  output_bam_file = output_dir + prefix + ".marked_duplicates.bam"
  output_metrics_file = output_bam_file + ".metrics.txt"
  command = ["gatk MarkDuplicates"]
  command += ["-I", bam_file]
  command += ["-O", output_bam_file]
  command += ["-M", output_metrics_file]
  command += ["--QUIET", "true"]
  command += ["--VERBOSITY", "ERROR"]
  command = " ".join(command)
  os.system(command)


def gatk_build_bam_index(io_tuple):

  bam_file, output_dir = io_tuple

  command = ["gatk BuildBamIndex"]
  command += ["-I", bam_file]
  command = " ".join(command)
  os.system(command)


def gatk_per_chromosome(chromosome):
  command_list = ["gatk --java-options '-Xmx2G'"]
  command_list += ["HaplotypeCaller"]
  # ~ command_list.append("Mutect2")
  # ~ command_list.append("--tumor-sample NIPT_600")
  command_list += ["--reference", hg38_fasta]
  command_list += ["--intervals", chromosome]
  command_list += ["--input", "temp_bam.list"]
  command_list += ["--output", "temp_" + chromosome + ".vcf.gz"]
  command = " ".join(command_list)
  os.system(command)


def gatk_multi_chromosome(bam_dir, output_file, num_threads=12):
  """Parallel gatk on multiple chromosomes via multiprocessing threads.

     Usage:
       bam_dir = "/data/nh2tran/GeneSolutions/temp/bam/"
       output_file = "/data/nh2tran/GeneSolutions/temp/sample.vcf.gz"
       gatk_multi_chromosome(bam_path, output_file, num_threads=8)
  """

  print("".join(["="] * 80)) # section-separating line
  print("gatk_multi_chromosome()")
  print("bam_dir =", bam_dir)
  print("output_file =", output_file)
  print("num_threads =", num_threads)
  print("".join(["="] * 80)) # section-separating line

  os.system("ls -1 {0}*.bam > temp_bam.list".format(bam_dir))
  chr_list = ["chr" + str(x) for x in range(1, 22+1)] + ["chrX", "chrY"]
  chr_batch = [chr_list[x:x+num_threads] for x in range(0, len(chr_list), num_threads)]
  for batch in chr_batch:
    pool = multiprocessing.Pool(processes=len(batch))
    pool.map(gatk_per_chromosome, batch)
    pool.close()
    pool.join()
  temp_vcf_list = ["temp_" + x + ".vcf.gz" for x in chr_list]
  os.system("bcftools concat {0} -Oz > {1}".format(" ".join(temp_vcf_list), output_file))
  os.system("rm temp_chr*")
  os.system("rm temp_bam*")


def extract_sample_variants_per_chromosome(argument):

  chromosome = argument[0]
  input_1000genomes = argument[1]
  population = argument[2]
  sample_file = argument[3]
  output_dir = argument[4]

  in_vcf = input_1000genomes + "ALL." + chromosome + "_GRCh38.genotypes.20170504.vcf.gz"
  out_vcf = output_dir + population + "." + chromosome + ".vcf.gz"
  # extract sample variants, update AC, AN
  command_subset = ["bcftools view"]
  command_subset += [in_vcf]
  command_subset += ["--samples-file", sample_file]
  command_subset += ["--min-ac=1"] # only variants that have non-reference allele in the subset samples
  #command_subset += ["--private"] # variants carried only by the subset samples
  #command_subset += ["--no-update"] # no update AC, AN
  command_subset += ["--force-samples"]
  command_subset = " ".join(command_subset)
  # update allele frequency
  command_af = "bcftools +fill-tags -- -t AF"
  # drop genotypes and output
  command_output = ["bcftools view"]
  command_output += ["--drop-genotypes"]
  command_output += ["-Oz >", out_vcf]
  command_output = " ".join(command_output)
  # stream 3 commands to 1 pipe
  command = " | ".join([command_subset, command_af, command_output])
  os.system(command)
  os.system("tabix -p vcf " + out_vcf)


def extract_sample_variants(input_1000genomes, population, sample_file, output_dir, num_threads=12):
  """Extract variants of population-specific samples from 1000 genomes project.
     Update AC, AN, AF.
     Drop sample genotypes, chromosome-merge to one vcf.

     Usage:
       input_1000genomes = "/data/nh2tran/GeneSolutions/1000genomes/release_20130502_GRCh38_positions/"
       population = "khv"
       sample_file = "/data/nh2tran/GeneSolutions/1000genomes/khv/khv_samples.txt"
       output_dir = "/data/nh2tran/GeneSolutions/1000genomes/temp/"
       extract_sample_variants(input_1000genomes, population, sample_file, output_dir)
  """
 
  print("".join(["="] * 80)) # section-separating line
  print("extract_sample_variants()")
  print("input_1000genomes =", input_1000genomes)
  print("population =", population)
  print("sample_file =", sample_file)
  print("output_dir =", output_dir)
  print("num_threads =", num_threads)
  print("".join(["="] * 80)) # section-separating line

  # group I/O arguments to batch to run parallel per chromosome
  chr_list = ["chr" + str(x) for x in range(1, 22+1)] + ["chrX", "chrY"]
  input_1000genomes_list = [input_1000genomes] * len(chr_list)
  population_list = [population] * len(chr_list)
  sample_file_list = [sample_file] * len(chr_list)
  output_dir_list = [output_dir] * len(chr_list)
  argument_list = list(zip(chr_list, input_1000genomes_list, population_list, sample_file_list, output_dir_list))
  argument_batch = [argument_list[x:x+num_threads] for x in range(0, len(argument_list), num_threads)]
  for batch in argument_batch:
    pool = multiprocessing.Pool(processes=len(batch))
    pool.map(extract_sample_variants_per_chromosome, batch)
    pool.close()
    pool.join()

  # chromosome-merge to one vcf
  vcf_list = [output_dir + population + "." + chromosome + ".vcf.gz" for chromosome in chr_list]
  out_vcf = output_dir + population + ".vcf.gz"
  command = ["bcftools concat"]
  command += vcf_list
  command += ["-Oz >", out_vcf]
  command = " ".join(command)
  os.system(command)
  os.system("tabix -p vcf " + out_vcf)


def read_vcf_stats(input_file):

  print("".join(["="] * 80)) # section-separating line
  print("read_vcf_stats()")
  print("".join(["="] * 80)) # section-separating line

  print("input_file:", input_file)
  af_list = []
  qual_list = []
  dp_list = []
  with open(input_file, 'r') as handle:
    for line in handle.readlines():
      if line[:2] == 'AF':
        af_list.append(line)
      elif line[:4] == 'QUAL':
        qual_list.append(line)
      elif line[:2] == 'DP':
        dp_list.append(line)
  print("len(af_list):", len(af_list))
  print("len(qual_list):", len(qual_list))
  print("len(dp_list):", len(dp_list))

  af_dup_dict = {}
  for af in af_list:
    af = af.strip().split('\t')
    af_id = int(af[1])
    af_value = int(float(af[2]) * 100)
    af_snp = int(af[3])
    af_dup_list = [af_value] * af_snp
    if af_id not in af_dup_dict:
      af_dup_dict[af_id] = af_dup_list
    else:
      af_dup_dict[af_id] += af_dup_list
  print("Number of SNPs:")
  for sample_id in af_dup_dict:
    print("Sample {0}: {1}".format(sample_id, len(af_dup_dict[sample_id])))
  print()
    
  qual_dup_dict = {}
  for qual in qual_list:
    qual = qual.strip().split('\t')
    qual_id = int(qual[1])
    qual_value = int(qual[2])
    qual_snp = int(qual[3])
    qual_indel = int(qual[6])
    qual_num = qual_snp + qual_indel
    qual_dup_list = [qual_value] * qual_num
    if qual_id not in qual_dup_dict:
      qual_dup_dict[qual_id] = qual_dup_list
    else:
      qual_dup_dict[qual_id] += qual_dup_list
  print("Number of SNPs and indels:")
  for sample_id in qual_dup_dict:
    print("Sample {0}: {1}".format(sample_id, len(qual_dup_dict[sample_id])))
  print()
    
  dp_dup_dict = {}
  for dp in dp_list:
    dp = dp.strip().split('\t')
    dp_id = int(dp[1])
    if dp[2] == '>500':
      dp[2] = '501'
    dp_value = int(dp[2])
    dp_num = int(dp[5])
    dp_dup_list = [dp_value] * dp_num
    if dp_id not in dp_dup_dict:
      dp_dup_dict[dp_id] = dp_dup_list
    else:
      dp_dup_dict[dp_id] += dp_dup_list
  print("Number of variant sites:")
  for sample_id in dp_dup_dict:
    print("Sample {0}: {1}".format(sample_id, len(dp_dup_dict[sample_id])))
  print()

  return af_dup_dict, qual_dup_dict, dp_dup_dict


# compare.sample_1_2.stats
# ~ af_dict, qual_dict, dp_dict = read_vcf_stats("/data/nh2tran/GeneSolutions/1000genomes/khv/khv.update_AC_AN_AF.vcf.gz.stats")
# ~ print("qual_dict[0] median:", np.median(qual_dict[0]))
# ~ print("qual_dict[2] median:", np.median(qual_dict[2]))
# ~ print("dp_dict[0] median:", np.median(dp_dict[0]))
# ~ print("dp_dict[2] median:", np.median(dp_dict[2]))
# ~ # boxplot
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([qual_dict[0], qual_dict[2]], labels=['Sample 1 only', 'Sample 1 and 2'])
# ~ ax.set_ylabel('Calling QUAL')
# ~ pyplot.savefig("compare.sample_1_2.qual.png")
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([dp_dict[0], dp_dict[2]], labels=['Sample 1 only', 'Sample 1 and 2'])
# ~ ax.set_ylim(0, 100)
# ~ ax.set_ylabel('DP')
# ~ pyplot.savefig("compare.sample_1_2.dp.png")

