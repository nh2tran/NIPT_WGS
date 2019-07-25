from __future__ import print_function


import os
import re
import math
import numpy as np
import pandas as pd
import csv
import multiprocessing
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
matplotlib.rcParams.update({'font.size': 11})

from Bio import SeqIO
from Bio.SeqIO import FastaIO


def get_AD_AF(isec_dir):
  """Get allele frequency from 1000 genomes or estimate from alt/ref read count.

     Usage:
       input_file =  "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/"
       get_AD_AF(isec_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("get_AD_AF()")
  print("isec_dir =", isec_dir)
  print("".join(["="] * 80)) # section-separating line

  vcf_1 = isec_dir + "0002.vcf"
  vcf_2 = isec_dir + "0003.vcf"
  with open(vcf_1, 'r') as handle_1:
    with open(vcf_2, 'r') as handle_2:
      # skip header lines
      # the number of header lines in two vcf files may be different
      for line_1 in handle_1:
        if line_1[0] != '#':
          break
      for line_2 in handle_2:
        if line_2[0] != '#':
          break
      
      # read a variant site and its pair of AD, AF
      num_sites = 0
      num_alleles = 0
      AD_freq_list = []
      AF_freq_list = []
      while line_1 and line_2:
        num_sites += 1
        chr_1, position_1, _, ref_1, alt_1, _, _, _, _, sample_1 = line_1.split('\t')
        chr_2, position_2, _, ref_2, alt_2, _, _, info_2 = line_2.split('\t')
        # make sure the variant site in two vcf files are consistent
        assert chr_1 == chr_2
        assert position_1 == position_2
        assert ref_1 == ref_2
        # order of alt alleles may be not consistent: "T,TA,TAAA" versus "TAAA,TA,T"
        alt_1 = alt_1.split(',')
        alt_2 = alt_2.split(',')
        assert set(alt_1) == set(alt_2)
        num_alleles += len(alt_1)

        
        # ref_count,alt_count[,alt_count]
        AD = sample_1.split(':')[1]
        AD_count = AD.split(',')
        AD_count = [int(x) for x in AD_count]
        AD_total_count = float(sum(AD_count))
        # exclude ref_count
        AD_count = AD_count[1:]
        # pair count to allele to avoid inconsistent order 
        AD_dict = dict(zip(alt_1, AD_count))
        
        # alt_freq[,alt_freq]
        AF = [x.split('=')[1] for x in info_2.split(';') if x[:3] == 'AF=']
        AF = AF[0]
        AF_freq = AF.split(',')
        AF_freq = [float(x) for x in AF_freq]
        # pair freq to allele to avoid inconsistent order 
        AF_dict = dict(zip(alt_2, AF_freq))
        
        # for each allele, AD_freq is estimated from alt/ref read count
        # AF_freq is from 1000 genomes
        assert len(AD_count) == len(AF_freq)
        for allele in AD_dict:
          AD_freq_list.append(AD_dict[allele] / AD_total_count)
          AF_freq_list.append(AF_dict[allele])
          
        # read next variant pair
        line_1 = handle_1.readline()
        line_2 = handle_2.readline()

  print("num_sites =", num_sites)
  print("num_alleles =", num_alleles)
  print("len(AD_freq_list) =", len(AD_freq_list))
  print("len(AF_freq_list) =", len(AF_freq_list))
  AD_file = isec_dir + "AD_freq_list.npy"
  AF_file = isec_dir + "AF_freq_list.npy"
  np.save(AD_file, AD_freq_list)
  np.save(AF_file, AF_freq_list)
  print(AD_file)
  print(AF_file)


def compare_AD_AF(AD_file, AF_file):
  """Compare allele frequency from 1000 genomes and estimation from alt/ref read count.

     Usage:
       AD_file =   "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02//AD_freq_list.npy"
       AF_file =   "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02//AF_freq_list.npy"
       compare_AD_AF(AD_file, AF_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("compare_AD_AF()")
  print("AD_file =", AD_file)
  print("AF_file =", AF_file)
  print("".join(["="] * 80)) # section-separating line

  AD_freq_list = np.load(AD_file)
  AF_freq_list = np.load(AF_file)
  print("AD_freq_list.shape =", AD_freq_list.shape)
  print("AF_freq_list.shape =", AF_freq_list.shape)
  corr = np.corrcoef(AD_freq_list, AF_freq_list)[0, 1]
  mse = np.mean((AD_freq_list - AF_freq_list) ** 2)
  median_ae = np.median(np.absolute(AD_freq_list - AF_freq_list))
  median_ape = np.median(np.absolute(AD_freq_list - AF_freq_list) / ((AD_freq_list + AF_freq_list) / 2))
  print("corr =", corr)
  print("mse_sqrt =", np.sqrt(mse))
  print("median_ae =", median_ae)
  print("median_ape =", median_ape)

  down_sampling = np.random.choice(AF_freq_list.shape[0], size=AF_freq_list.shape[0]//100, replace=False)
  fig, ax = pyplot.subplots(figsize=(8, 8))
  pyplot.scatter(x=AF_freq_list[down_sampling], y=AD_freq_list[down_sampling],
                 marker='.', s=0.1)
  ax.set_xlabel('KHV AF frequency')
  ax.set_ylabel('gatk AD frequency')
  pyplot.show()


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


def mapping_bwa(index, fastq_dir, output_dir, read_group_dict, num_threads=1, num_test=1):
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


def process_multi_samples(function, bam_list_file, output_dir, num_threads=1, num_test=1):
  """Process multiple samples in parallel threads.

     Usage:
       function = samtools_flagstat_per_bam
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam_list.txt"
       output_dir = "/data/nh2tran/GeneSolutions/temp/samtools_flagstat/"
       num_threads = 1
       num_test = 1
       process_multi_samples(function, bam_list_file, output_dir, num_threads, num_test)
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


def samtools_flagstat(bam_list_file, output_dir, num_threads=1, num_test=1):
  """Summarize bam files.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam_list.txt"
       output_dir = "/data/nh2tran/GeneSolutions/temp/samtools_flagstat/"
       samtools_flagstat(bam_list_file, output_dir, num_threads=1, num_test=1)
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


def filter_2plus_reads(io_tuple):

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]
  prefix = bam_file_name[:-4]

  bed_file = output_dir + bam_file_name + ".2plus.bed"
  command = ["bedtools genomecov"]
  command += ["-bg"]
  command += ["-ibam", bam_file]
  command += ["| grep -vw \"1$\""]
  command += [">", bed_file]
  command = " ".join(command)
  os.system(command)

  filtered_bam_file = output_dir + prefix + ".filtered_2plus.bam"
  command = ["bedtools intersect"]
  command += ["-v"]
  command += ["-abam", bam_file]
  command += ["-b", bed_file]
  command += [">", filtered_bam_file]
  command = " ".join(command)
  os.system(command)

  command = ["gatk BuildBamIndex"]
  command += ["-I", filtered_bam_file]
  command = " ".join(command)
  os.system(command)


def genomecov_hist_per_bam(io_tuple):

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]

  hist_file = output_dir + bam_file_name + ".hist.genome"
  command = ["bedtools genomecov"]
  command += ["-ibam", bam_file]
  command += ["| grep \"genome\""]
  command += [">", hist_file]
  command = " ".join(command)
  os.system(command)


def genomecov_hist(bam_list_file, output_dir, num_threads=1, num_test=1):
  """Calculate bedtools genomecov on a list of bam files.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/genomecov_hist/"
       genomecov_hist(bam_list_file, output_dir, num_threads=1, num_test=1)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_hist()")
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("".join(["="] * 80)) # section-separating line

  # run bedtools genomecov to get hist and bed files
  # took ~5 hours
  #process_multi_samples(genomecov_hist_per_bam, bam_list_file, output_dir, num_threads=num_threads, num_test=num_test)

  hist_file_list = [output_dir + x for x in os.listdir(output_dir) if "hist.genome" in x]
  print("len(hist_file_list) =", len(hist_file_list))

  frequency_zero_4plus= []
  mean_depth = []
  for hist_file in hist_file_list:
    hist_array = np.loadtxt(hist_file, usecols=(1,2))
    depth = hist_array[:,0]
    frequency = hist_array[:,1]
    mean_depth.append(np.average(depth, weights=frequency))
    # record frequency of zero, 1x, 2x, 3x, 4plus coverage
    max_depth = len(frequency) - 1
    if max_depth < 4:
      frequency = np.append(frequency, [0]*(4-max_depth))
    frequency[4] = np.sum(frequency[4:])
    frequency_zero_4plus.append(frequency[:5])
  frequency_zero_4plus = np.array(frequency_zero_4plus)
  mean_depth = np.array(mean_depth)

  ref_length = frequency_zero_4plus[0,:].sum()
  print("ref_length =", ref_length)
  # genome coverage per bam
  pct_nonzero = np.mean((ref_length - frequency_zero_4plus[:,0]) / ref_length)
  pct_1x = np.mean(frequency_zero_4plus[:,1] / ref_length)
  pct_2x = np.mean(frequency_zero_4plus[:,2] / ref_length)
  pct_3x = np.mean(frequency_zero_4plus[:,3] / ref_length)
  pct_4plus = np.mean(frequency_zero_4plus[:,4] / ref_length)
  print("pct_nonzero =", pct_nonzero)
  print("pct_1x =", pct_1x)
  print("pct_2x =", pct_2x)
  print("pct_3x =", pct_3x)
  print("pct_4plus =", pct_4plus)
  # depth per bam
  print("mean_depth =", np.average(mean_depth))

  # calculate #bam at depth-2plus genome positions 
  # ~ bed_file_list = [output_dir + x for x in os.listdir(output_dir) if "2plus.bed" in x]
  # ~ print("len(bed_file_list) =", len(bed_file_list))
  # use {chr_name: array} to represent genome positions 
  # ~ genome_array = {}
  # ~ chr_list = ['chr' + str(x) for x in range(1, 23)]
  # ~ chr_list += ['chrX', 'chrY']
  # ~ with open("hg38_selected.genome", 'r') as file_handle:
    # ~ for line in file_handle:
      # ~ chr_name, chr_length = line.split()
      # ~ chr_length = int(chr_length)
      # ~ genome_array[chr_name]  = np.zeros(chr_length, dtype='int16')
  # ~ ref_length = sum([len(x) for x in genome_array.values()])
  # ~ print("ref_length =", ref_length)
  # ~ for bed_file in bed_file_list:
    # ~ print(bed_file)
    # ~ with open(bed_file, 'r') as file_handle:
      # ~ for line in file_handle:
        # ~ chr_name, start, end, depth = line.split()
        # ~ start, end, depth = [int(x) for x in [start, end, depth]]
        # ~ if depth == 2:
          # ~ genome_array[chr_name][start:end] += 1
  # ~ os.system("mkdir genome_array")
  # ~ for chr_name in genome_array:
    # ~ np.save("genome_array/" + chr_name, genome_array[chr_name])
  # ~ genome_array_dir = "/data/nh2tran/GeneSolutions/temp/genome_array/"
  # ~ genome_array = {}
  # ~ chr_list = ['chr' + str(x) for x in range(1, 23)]
  # ~ chr_list += ['chrX', 'chrY']
  # ~ for chr_name in chr_list[:1]:
      # ~ genome_array[chr_name] = np.load(genome_array_dir + chr_name + ".npy")
      # ~ print(chr_name)
      # ~ print(len(genome_array[chr_name]))
      # ~ print(np.sum(genome_array[chr_name] >= 1))
      # ~ print(np.max(genome_array[chr_name]))
      # ~ print(np.median(genome_array[chr_name]))
      # ~ #pyplot.hist(genome_array[chr_name], bins=range(1, 51))
      # ~ #pyplot.show()
      # ~ with open("chr1_depth_2x_samples.bed", 'w') as file_handle:
          # ~ for index in np.flatnonzero(genome_array[chr_name]):
              # ~ line = '\t'.join(['chr1', str(index), str(index+1), str(genome_array[chr_name][index])])
              # ~ line += '\n'
              # ~ file_handle.write(line)


def genomecov_merge(hist_file):
  """Calculate bedtools genomecov on a merged bam file.

     Usage:
       hist_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.hist.genome"
       genomecov_merge(hist_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_merge()")
  print("hist_file =", hist_file)
  print("".join(["="] * 80)) # section-separating line
  
  hist_array = np.loadtxt(hist_file, usecols=(1,2))
  depth = hist_array[:,0]
  frequency = hist_array[:,1]
  
  ref_length = np.sum(frequency)
  print("ref_length =", ref_length)
  
  # depth mean/median
  mean_depth = np.average(depth, weights=frequency)
  print("mean_depth =", mean_depth)
  frequency_cumsum = np.cumsum(frequency)
  median_index = np.searchsorted(frequency_cumsum, frequency_cumsum[-1]/2.0)
  median_depth = depth[median_index]
  print("median_depth =", median_depth)
  
  # depth sd
  error2 = np.square(depth - mean_depth)
  mean_sd = np.sqrt(np.average(error2, weights=frequency))
  print("mean_sd =", mean_sd)
  sorted_indices = np.argsort(error2)
  error2_sorted = error2[sorted_indices]
  frequency_sorted = frequency[sorted_indices]
  frequency_sorted_cumsum = np.cumsum(frequency_sorted)
  median_index = np.searchsorted(frequency_sorted_cumsum, frequency_sorted_cumsum[-1]/2.0)
  median_sd = np.sqrt(error2_sorted[median_index])
  print("median_sd =", median_sd)
  
  # depth percentiles
  depth_zero_pct = frequency[0] / frequency_cumsum[-1]
  print("depth_zero_pct =", depth_zero_pct)
  print("depth_nonzero_pct =", 1-depth_zero_pct)
  pct_10_depth = np.flatnonzero(frequency_cumsum / frequency_cumsum[-1] >= 0.1)[0]
  pct_99_depth = np.flatnonzero(frequency_cumsum / frequency_cumsum[-1] >= 0.99)[0]
  print("pct_10_depth =", pct_10_depth)
  print("pct_99_depth =", pct_99_depth)
  print("max_depth", depth[-1])
  

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


def create_chromosome_group(ref_dict_file):

  # This is adopted from gatk scatter.
  # Basically chromosomes are grouped so that the group total lengths are similar for efficient parallelization.
  # Create list of sequences for scatter-gather parallelization
  # call CreateSequenceGroupingTSV{}
  with open(ref_dict_file, "r") as file_handle:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in file_handle:
        if line.startswith("@SQ"):
          line_split = line.split("\t")
          # (Sequence_Name, Sequence_Length)
          sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
  # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
  # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
  hg38_protection_tag = ":1+"
  # initialize the tsv string with the first sequence
  tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
  temp_size = sequence_tuple_list[0][1]
  for sequence_tuple in sequence_tuple_list[1:]:
      if temp_size + sequence_tuple[1] <= longest_sequence:
          temp_size += sequence_tuple[1]
          tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
      else:
          tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
          temp_size = sequence_tuple[1]
  # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
  tsv_string += '\n' + "unmapped"

  # return list of groups instead of writing to files as in gatk wdl pipeline
  sequence_grouping_with_unmapped = [group.split('\t') for group in tsv_string.split('\n')]
  sequence_grouping = sequence_grouping_with_unmapped[:-1]
  return sequence_grouping, sequence_grouping_with_unmapped


def gatk_BaseRecalibrator_per_intervals(args):

  bam_file, ref_file, dbSNP_vcf, intervals, output_file = args
  # Generate the recalibration model by interval
  command = ["gatk BaseRecalibrator"]
  command += ["-I", bam_file]
  command += ["-R", ref_file]
  command += ["--known-sites", dbSNP_vcf]
  command += ["-L", " -L ".join(intervals)]
  command += ["-O", output_file]
  command += ["--use-original-qualities"]
  command = " ".join(command)
  os.system(command)


def gatk_ApplyBQSR_per_intervals(args):

  bam_file, ref_file, intervals, bqsr_report_file, output_file = args
  # Generate the recalibration model by interval
  command = ["gatk ApplyBQSR"]
  command += ["-I", bam_file]
  command += ["-R", ref_file]
  command += ["-L", " -L ".join(intervals)]
  command += ["-bqsr", bqsr_report_file]
  command += ["-O", output_file]
  command += ["--use-original-qualities"]
  command = " ".join(command)
  os.system(command)


def gatk_bqsr_workflow(bam_list_file, ref_dict_file, ref_file, dbSNP_vcf, output_dir, num_threads, num_test=1):
  """Usage
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       ref_dict_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.dict"
       ref_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
       dbSNP_vcf = "/data/nh2tran/GeneSolutions/dbSNP/human_9606_b150_GRCh38p7/All_20170710.with_chr.vcf.gz"
       output_dir = "/data/nh2tran/GeneSolutions/temp/gatk_bqsr_workflow/"
       num_threads = 16
       gatk_bqsr_workflow(bam_file_list, ref_dict_file, ref_file, dbSNP_vcf, output_dir, num_threads)
  """

  # create a folder to store results
  if output_dir is not None:
    assert not os.path.exists(output_dir), "output_dir exists " + output_dir
    os.makedirs(output_dir)

  # Create groups of chromosomes for parallelization.
  # This is adopted from gatk scatter.
  # Basically chromosomes are grouped so that the group total lengths are similar for efficient parallelization.
  chr_group_pool, chr_group_with_unmapped_pool = create_chromosome_group(ref_dict_file)

  with open(bam_list_file, 'r') as file_handle:
    bam_list = [line.strip() for line in file_handle.readlines()]
  print("Number of bam files:", len(bam_list))
  # limite the number of bam files in test mode
  if num_test is not None:
    bam_list = bam_list[:num_test]

  for bam_file in bam_list:

    bam_basename = bam_file.split('/')[-1][:-4] # remove path and .bam from file name
    bqsr_basename = output_dir + "bqsr" + "." + bam_basename

    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel chromosome groups
    bam_file_pool = [bam_file] * len(chr_group_pool)
    ref_file_pool = [ref_file] * len(chr_group_pool)
    dbSNP_vcf_pool = [dbSNP_vcf] * len(chr_group_pool)
    bqsr_report_pool = [bqsr_basename + "." + "_".join(chr_group) + ".recal_data.csv" for chr_group in chr_group_pool]
    args_pool = zip(bam_file_pool, ref_file_pool, dbSNP_vcf_pool, chr_group_pool, bqsr_report_pool)
    pool = multiprocessing.Pool(processes=num_threads)
    pool.map(gatk_BaseRecalibrator_per_intervals, args_pool)
    pool.close()
    pool.join()


    # Merge the recalibration reports resulting from by-interval recalibration
    bqsr_report_file = bqsr_basename + ".recal_data.csv"
    command = ["gatk GatherBQSRReports"]
    command += ["-I", " -I ".join(bqsr_report_pool)]
    command += ["-O", bqsr_report_file]
    command = " ".join(command)
    os.system(command)

    # Apply the recalibration model by interval
    bam_file_pool = [bam_file] * len(chr_group_with_unmapped_pool)
    ref_file_pool = [ref_file] * len(chr_group_with_unmapped_pool)
    bqsr_report_file_pool = [bqsr_report_file] * len(chr_group_with_unmapped_pool)
    bam_pool = [bqsr_basename + "." + "_".join(chr_group) + ".recalibrated.bam" for chr_group in chr_group_with_unmapped_pool]
    args_pool = zip(bam_file_pool, ref_file_pool, chr_group_with_unmapped_pool, bqsr_report_file_pool, bam_pool)
    pool = multiprocessing.Pool(processes=num_threads)
    pool.map(gatk_ApplyBQSR_per_intervals, args_pool)
    pool.close()
    pool.join()

    # Merge the recalibrated BAM files resulting from by-interval recalibration
    output_bam_file = output_dir + bam_basename + ".bqsr.bam"
    command = ["gatk GatherBamFiles"]
    command += ["-I", " -I ".join(bam_pool)]
    command += ["-O", output_bam_file]
    command += ["--CREATE_INDEX", "true"]
    command = " ".join(command)
    os.system(command)

    # Remove temporary bsqr files
    os.system("rm " + output_dir + "bqsr.*")


def gatk_per_chromosome(args):

  bam_list_file, reference_file, chromosome = args

  command_list = ["gatk --java-options '-Xmx2G'"]

  # ~ command_list += ["HaplotypeCaller"]
  # Ploidy (number of chromosomes) per sample.
  # For pooled data, set to (Number of samples in each pool * Sample Ploidy).
  # Default value: 2. 
  # ~ command_list += ["--sample-ploidy", str(10)]

  command_list.append("Mutect2")
  command_list.append("--tumor-sample NIPT_600")

  command_list += ["--input", bam_list_file]
  command_list += ["--reference", reference_file]
  command_list += ["--intervals", chromosome]
  command_list += ["--output", "temp_" + chromosome + ".vcf.gz"]
  command = " ".join(command_list)
  os.system(command)


def gatk_multi_chromosome(bam_list_file, reference_file, output_file, num_threads=1):
  """Parallel gatk on multiple chromosomes via multiprocessing threads.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       reference_file = "/data/nh2tran/GeneSolutions/references/hg38_selected.fa"
       output_file = "/data/nh2tran/GeneSolutions/temp/sample.vcf.gz"
       gatk_multi_chromosome(bam_path, output_file, num_threads=8)
  """

  print("".join(["="] * 80)) # section-separating line
  print("gatk_multi_chromosome()")
  print("bam_list_file =", bam_list_file)
  print("output_file =", output_file)
  print("num_threads =", num_threads)
  print("".join(["="] * 80)) # section-separating line

  chr_list = ["chr" + str(x) for x in range(1, 22+1)] + ["chrX", "chrY"]
  bam_list = [bam_list_file] * len(chr_list)
  reference_list = [reference_file] * len(chr_list)
  args_list = zip(bam_list, reference_list, chr_list)
  pool = multiprocessing.Pool(processes=num_threads)
  pool.map(gatk_per_chromosome, args_list)
  pool.close()
  pool.join()
  temp_vcf_list = ["temp_" + x + ".vcf.gz" for x in chr_list]
  os.system("bcftools concat {0} -Oz > {1}".format(" ".join(temp_vcf_list), output_file))
  os.system("rm temp_chr*")


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
  command_subset += ["--min-ac 1"] # only variants that have non-reference allele in the subset samples
  command_subset += ["--trim-alt-alleles"] # from multi-allele sites: C,T AC=198,0;AF=1,0;
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


def extract_sample_variants(input_1000genomes, population, sample_file, output_dir, num_threads):
  """Extract variants of population-specific samples from 1000 genomes project.
     Update AC, AN, AF.
     Drop sample genotypes, chromosome-merge to one vcf.

     Usage:
       input_1000genomes = "/data/nh2tran/GeneSolutions/1000genomes/release_20130502_GRCh38_positions/"
       population = "khv"
       sample_file = "/data/nh2tran/GeneSolutions/1000genomes/khv/khv_samples.txt"
       output_dir = "/data/nh2tran/GeneSolutions/1000genomes/temp/"
       num_threads = 1
       extract_sample_variants(input_1000genomes, population, sample_file, output_dir, num_threads)
  """
 
  print("".join(["="] * 80)) # section-separating line
  print("extract_sample_variants()")
  print("input_1000genomes =", input_1000genomes)
  print("population =", population)
  print("sample_file =", sample_file)
  print("output_dir =", output_dir)
  print("num_threads =", num_threads)
  print("".join(["="] * 80)) # section-separating line

  # run parallel per chromosome
  chr_list = ["chr" + str(x) for x in range(1, 22+1)] + ["chrX", "chrY"]
  input_1000genomes_list = [input_1000genomes] * len(chr_list)
  population_list = [population] * len(chr_list)
  sample_file_list = [sample_file] * len(chr_list)
  output_dir_list = [output_dir] * len(chr_list)
  argument_list = list(zip(chr_list, input_1000genomes_list, population_list, sample_file_list, output_dir_list))
  pool = multiprocessing.Pool(processes=num_threads)
  pool.map(extract_sample_variants_per_chromosome, argument_list)
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
# ~ vcf_file = "/data/nh2tran/GeneSolutions/temp/compare.NIPT_100.HaplotypeCaller.khv.stats"
# ~ af_dict, qual_dict, dp_dict = read_vcf_stats(vcf_file)
# ~ print("qual_dict[0] median:", np.median(qual_dict[0]))
# ~ print("qual_dict[2] median:", np.median(qual_dict[2]))
# ~ print("dp_dict[0] median:", np.median(dp_dict[0]))
# ~ print("dp_dict[2] median:", np.median(dp_dict[2]))
# ~ # boxplot
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([qual_dict[0], qual_dict[2]], labels=['Sample 0 only', 'Sample 0 and 2'])
# ~ ax.set_ylabel('Calling QUAL')
# ~ pyplot.show()
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([dp_dict[0], dp_dict[2]], labels=['Sample 0 only', 'Sample 0 and 2'])
# ~ ax.set_ylim(0, 100)
# ~ ax.set_ylabel('DP')


# ~ af_dict, qual_dict, dp_dict = read_vcf_stats("/data/nh2tran/GeneSolutions/1000genomes/khv/khv.update_AC_AN_AF.vcf.gz.stats")
# ~ # AF summary
# ~ print("AF mean:", np.mean(af_dict[0]))
# ~ print("AF median:", np.median(af_dict[0]))
# ~ # boxplot
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([af_dict[0]], labels=['Sample 0'])
# ~ ax.set_ylabel('Allele frequency')
# ~ pyplot.savefig("af_boxplot.png")
# ~ # histogram
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.hist([af_dict[0]], bins=20, weights=np.ones_like(af_dict[0]) / float(len(af_dict[0])))
# ~ ax.set_xlabel('Allele frequency')
# ~ ax.set_ylabel('Percentage of SNPs')
# ~ pyplot.savefig("af_hist.png")


# ~ # QUAL summary
# ~ print("QUAL mean:", np.mean(qual_dict[0]))
# ~ print("QUAL median:", np.median(qual_dict[0]))
# ~ # boxplot
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([qual_dict[0]], labels=['Sample 0'])
# ~ ax.set_ylabel('Calling QUAL')
# ~ pyplot.savefig("qual_boxplot.png")
# ~ # histogram
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.hist([qual_dict[0]], bins=20, weights=np.ones_like(qual_dict[0]) / float(len(qual_dict[0])))
# ~ ax.set_xlabel('Calling QUAL')
# ~ ax.set_ylabel('Percentage of variants')
# ~ pyplot.savefig("qual_hist.png")


# ~ # DP summary
# ~ print("DP mean:", np.mean(dp_dict[0]))
# ~ print("DP median:", np.median(dp_dict[0]))
# ~ # boxplot
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.boxplot([dp_dict[0]], labels=['Sample 0'])
# ~ ax.set_ylabel('Calling DP')
# ~ pyplot.savefig("dp_boxplot.png")
# ~ # histogram
# ~ fig, ax = pyplot.subplots()
# ~ pyplot.hist([dp_dict[0]], bins=20, weights=np.ones_like(dp_dict[0]) / float(len(dp_dict[0])))
# ~ ax.set_xlabel('Calling DP')
# ~ ax.set_ylabel('Percentage of variants')
# ~ pyplot.savefig("dp_hist.png")

