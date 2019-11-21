from __future__ import print_function


import os
import sys
import re
import math
import numpy as np
import pandas as pd
import csv
import multiprocessing
import time
import pickle

import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
matplotlib.rcParams.update({'font.size': 11})

from Bio import SeqIO
from Bio.SeqIO import FastaIO


def preprocessing_per_sample(io_list):
  """fastqc | trim | align | sort | dedup | index
  """

  r1_r2, index, dir_tuple, read_group_dict = io_list
  r1, r2 = r1_r2
  fastqc_dir, trim_paired_dir, trim_unpaired_dir, bam_dir = dir_tuple
  sample_id = r1.split('/')[-1].split('_')[0]

  num_threads = 4
  java_xmx = "-Xmx16G"

  # fastqc
  command_fastqc = ["fastqc"]
  command_fastqc += [r1, r2]
  command_fastqc += ["-o", fastqc_dir]
  command_fastqc += ["-t", str(2)] # cannot be more than number of input files
  command_fastqc += ["--quiet"]
  command_fastqc = " ".join(command_fastqc)
  os.system(command_fastqc)

  # trim
  trim_r1_paired = trim_paired_dir + sample_id + ".R1.paired.fastq.gz"
  trim_r1_unpaired = trim_unpaired_dir + sample_id + ".R1.unpaired.fastq.gz"
  trim_r2_paired = trim_paired_dir + sample_id + ".R2.paired.fastq.gz"
  trim_r2_unpaired = trim_unpaired_dir + sample_id + ".R2.unpaired.fastq.gz"
  command_trim = ["java {0:s} -jar $trimmomatic PE".format(java_xmx)]
  command_trim += [r1, r2]
  command_trim += [trim_r1_paired, trim_r1_unpaired, trim_r2_paired, trim_r2_unpaired]
  command_trim += ["ILLUMINACLIP:/data/nh2tran/GeneSolutions/tools/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10"]
  command_trim += ["LEADING:3"]
  command_trim += ["TRAILING:3"]
  command_trim += ["SLIDINGWINDOW:4:15"]
  command_trim += ["CROP:75"]
  command_trim += ["MINLEN:36"]
  command_trim += ["-threads", str(num_threads)]
  command_trim = " ".join(command_trim)
  os.system(command_trim)

  # piping align | sort
  temp_file = bam_dir + sample_id + ".sorted.sam"
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
  command_bwa += ["-t", str(num_threads)]
  command_bwa = " ".join(command_bwa)
  command_sort = ["samtools sort"]
  command_sort += ["-O", "SAM"]
  command_sort += ["-@", str(num_threads)]
  command_sort += [">", temp_file]
  command_sort = " ".join(command_sort)
  pipe_bwa_sort = " | ".join([command_bwa, command_sort])
  os.system(pipe_bwa_sort)

  # dedup and index
  bam_file = bam_dir + sample_id + ".sorted.deduped.bam"
  metrics_file = bam_file + ".metrics.txt"
  command_dedup = ["gatk --java-options '{0:s}' MarkDuplicates".format(java_xmx)]
  command_dedup += ["-I", temp_file]
  command_dedup += ["-O", bam_file]
  command_dedup += ["-M", metrics_file]
  command_dedup += ["--CREATE_INDEX", "true"]
  command_dedup += ["--QUIET", "true"]
  command_dedup += ["--VERBOSITY", "ERROR"]
  command_dedup = " ".join(command_dedup)
  os.system(command_dedup)

  # remove temp files
  os.system("rm " + trim_r1_paired)
  os.system("rm " + trim_r2_paired)
  os.system("rm " + temp_file)


def preprocessing(index, fastq_dir, output_dir, read_group_dict, num_samples, num_test):
  """fastqc | trim | align | sort | dedup | index
     Parallel on num_samples

     Usage:
       index = "/data/nh2tran/GeneSolutions/indexes/hg38_selected"
       fastq_dir = "/data/nh2tran/GeneSolutions/fastq/NIPT_700/"
       output_dir = "/data/nh2tran/GeneSolutions/temp/preprocessing/"
       read_group_dict = {'SM': "NIPT_700", 'PL': "unknown", 'PU': "unknown"}
       num_samples = 1 # default 4 threads, 16G per sample
       num_test = 1
       preprocessing(index, fastq_dir, output_dir, read_group_dict, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("preprocessing()")
  print("index =", index)
  print("fastq_dir =", fastq_dir)
  print("output_dir =", output_dir)
  print("read_group_dict =", read_group_dict)
  print("num_samples =", num_samples)
  print("num_test =", num_test)
  print("".join(["="] * 80)) # section-separating line

  # get paired-end fastq and check if they match
  fastq_list = [os.path.join(fastq_dir, x) for x in os.listdir(fastq_dir)]
  r1_list = [x for x in fastq_list if 'R1' in x]
  r2_list = [x for x in fastq_list if 'R2' in x]
  r1_list.sort()
  r2_list.sort()
  if num_test is not None:
    r1_list = r1_list[:num_test]
    r2_list = r2_list[:num_test]
  for r1, r2 in zip(r1_list, r2_list):
    assert r1 == r2.replace('_R2_', '_R1_'), "Error: R1 and R2 not matched"
  print("Number of paired-end fastq:", len(r1_list))
  print()

  # prepare the output folders
  fastqc_dir = output_dir + "fastqc/"
  # after trimming, paired reads are kept separately to be removed after alignment
  trim_paired_dir = output_dir + "trim_paired/"
  trim_unpaired_dir = output_dir + "trim_unpaired/"
  bam_dir = output_dir + "bam/"
  dir_tuple = (fastqc_dir, trim_paired_dir, trim_unpaired_dir, bam_dir)
  os.system("mkdir " + output_dir)
  os.system("mkdir " + fastqc_dir)
  os.system("mkdir " + trim_paired_dir)
  os.system("mkdir " + trim_unpaired_dir)
  os.system("mkdir " + bam_dir)

  # preprocessing num_samples in parallel
  r1_r2_list = zip(r1_list, r2_list)
  args_list = [[r1_r2, index, dir_tuple, read_group_dict] for r1_r2 in r1_r2_list]
  pool = multiprocessing.Pool(processes=num_samples)
  pool.map(preprocessing_per_sample, args_list)
  pool.close()
  pool.join()


def process_multi_samples(function, bam_list_file, output_dir, num_samples, num_test):
  """Process multiple samples in parallel.

     Usage:
       function = samtools_flagstat_per_bam
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam_list.txt"
       output_dir = "/data/nh2tran/GeneSolutions/temp/samtools_flagstat/"
       num_samples = 1
       num_test = 1
       process_multi_samples(function, bam_list_file, output_dir, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("process_multi_samples()")
  print("function =", function.__name__)
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("num_samples =", num_samples)
  print("num_test =", num_test)
  print("".join(["="] * 80)) # section-separating line

  with open(bam_list_file, 'r') as file_handle:
    bam_list = [line.strip() for line in file_handle.readlines()]
  print("Number of bam files:", len(bam_list))
  # limite the number of bam files in test mode
  if num_test is not None:
    bam_list = bam_list[:num_test]

  output_dir_list = [output_dir] * len(bam_list)
  argument_list = zip(bam_list, output_dir_list)
  pool = multiprocessing.Pool(processes=num_samples)
  pool.map(function, argument_list)
  pool.close()
  pool.join()

  print("multiprocessing.Pool() finished.")


def genomecov_hist_per_bam(io_tuple):
  """Run bedtools genomecov hist on 1 bam file, and extract genome lines.
     bedtools can only use 1 thread.
  """

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]

  hist_file = output_dir + bam_file_name + ".genomecov_hist"
  command = ["bedtools genomecov"]
  command += ["-ibam", bam_file]
  command += ["| grep -w \"^genome\""]
  command += [">", hist_file]
  command = " ".join(command)
  os.system(command)


def genomecov_hist(bam_list_file, output_dir, num_samples, num_test):
  """Run bedtools genomecov hist on a list of bam files.
     Summarize percentage of nonzero, 1x, 2x, etc, and depth per bam.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/genomecov_hist/"
       num_samples = 1
       num_test = 1
       genomecov_hist(bam_list_file, output_dir, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_hist()")
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("num_samples =", num_samples)
  print("num_test =", num_test)
  print("".join(["="] * 80)) # section-separating line

  # run bedtools genomecov hist on num_samples in parallel
  os.system("mkdir " + output_dir)
  process_multi_samples(genomecov_hist_per_bam, bam_list_file, output_dir, num_samples, num_test)

  # summarize coverage from genomecov_hist
  hist_file_list = [output_dir + x for x in os.listdir(output_dir) if ".genomecov_hist" in x]
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


def genomecov_hist_sum(hist_file):
  """Summarize coverage and depth from a bedtools genomecov hist file.

     Usage:
       hist_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.genomecov_hist"
       genomecov_hist_sum(hist_file)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_hist_sum()")
  print("hist_file =", hist_file)
  print("".join(["="] * 80)) # section-separating line
  
  depth_dict = {}
  with open(hist_file, 'r') as filehandle:
    for line in filehandle:
      chr_name, depth, frequency = line.split()
      depth = int(depth)
      frequency = int(frequency)
      if chr_name in depth_dict:
        depth_dict[chr_name].append([depth, frequency])
      else:
        depth_dict[chr_name] = [[depth, frequency]]

  # genome depth
  depth = np.array([x[0] for x in depth_dict["genome"]])
  frequency = np.array([x[1] for x in depth_dict["genome"]])
  # depth mean/median
  mean_depth = np.average(depth, weights=frequency)
  frequency_cumsum = np.cumsum(frequency)
  median_index = np.searchsorted(frequency_cumsum, frequency_cumsum[-1]/2.0)
  median_depth = depth[median_index]
  # depth sd mean/median
  error2 = np.square(depth - mean_depth)
  mean_sd = np.sqrt(np.average(error2, weights=frequency))
  #error2 = np.square(depth - median_depth)
  sorted_indices = np.argsort(error2)
  error2_sorted = error2[sorted_indices]
  frequency_sorted = frequency[sorted_indices]
  frequency_sorted_cumsum = np.cumsum(frequency_sorted)
  median_index = np.searchsorted(frequency_sorted_cumsum, frequency_sorted_cumsum[-1]/2.0)
  median_sd = np.sqrt(error2_sorted[median_index])
  # depth percentiles
  depth_zero_pct = frequency[0] / frequency_cumsum[-1]
  pct_10_depth = np.flatnonzero(frequency_cumsum / frequency_cumsum[-1] >= 0.1)[0]
  pct_99_depth = np.flatnonzero(frequency_cumsum / frequency_cumsum[-1] >= 0.99)[0]

  print("mean_depth =", mean_depth)
  print("median_depth =", median_depth)
  print("mean_sd =", mean_sd)
  print("median_sd =", median_sd)
  print("depth_zero_pct =", depth_zero_pct)
  print("depth_nonzero_pct =", 1-depth_zero_pct)
  print("pct_10_depth =", pct_10_depth)
  print("pct_99_depth =", pct_99_depth)
  print("max_depth", depth[-1])

  depth_lower_bound = median_depth - 3 * median_sd
  depth_upper_bound = median_depth + 3 * median_sd #float('inf')
  depth_bound_count = sum([y for x, y in zip(depth, frequency) if (x >= depth_lower_bound and x <= depth_upper_bound)])
  depth_bound_pct = depth_bound_count / frequency_cumsum[-1]
  print("depth_lower_bound =", depth_lower_bound)
  print("depth_upper_bound =", depth_upper_bound)
  print("depth_bound_pct =", depth_bound_pct)

  # depth distribution
  png_file =  "plot_depth_distr.png"
  fig, ax = pyplot.subplots()
  pyplot.hist(depth, weights=frequency, bins=50, range=(0, 800),
              facecolor='salmon', alpha=0.8, rwidth=0.8)
  pyplot.ylabel("Number of genome positions")
  pyplot.xlabel("Sequencing depth")
  pyplot.text(475, 2.5e8, 'average=364, error=67')
  pyplot.grid(axis='y', alpha=0.8)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  pyplot.savefig(png_file)
  print(png_file)

  chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
  for chr_name in chr_list:
    depth = np.array([x[0] for x in depth_dict[chr_name]])
    frequency = np.array([x[1] for x in depth_dict[chr_name]])
    # depth mean/median
    mean_depth = np.average(depth, weights=frequency)
    frequency_cumsum = np.cumsum(frequency)
    median_index = np.searchsorted(frequency_cumsum, frequency_cumsum[-1]/2.0)
    median_depth = depth[median_index]
    # depth sd mean/median
    error2 = np.square(depth - mean_depth)
    mean_sd = np.sqrt(np.average(error2, weights=frequency))
    error2 = np.square(depth - median_depth)
    sorted_indices = np.argsort(error2)
    error2_sorted = error2[sorted_indices]
    frequency_sorted = frequency[sorted_indices]
    frequency_sorted_cumsum = np.cumsum(frequency_sorted)
    median_index = np.searchsorted(frequency_sorted_cumsum, frequency_sorted_cumsum[-1]/2.0)
    median_sd = np.sqrt(error2_sorted[median_index])
    print(chr_name, mean_depth, mean_sd, median_depth, median_sd)
  

def genomecov_2x_per_bam(io_tuple):
  """Run bedtools genomecov bg on 1 bam file, and extract 2x regions.
     bedtools can only use 1 thread.
  """

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]

  bed_file = output_dir + bam_file_name + ".genomecov_2x"
  command = ["bedtools genomecov"]
  command += ["-bg"]
  command += ["-ibam", bam_file]
  command += ["| grep -w \"2$\""]
  command += [">", bed_file]
  command = " ".join(command)
  os.system(command)


def genomecov_2x(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test):
  """Run bedtools genomecov bg on a list of bam files.
     Summarize distribution of 2x regions across genome and samples.

     Usage:
       hg38_selected_genome = "/data/nh2tran/GeneSolutions/references/hg38_selected.genome"
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/genomecov_2x/"
       num_samples = 1
       num_test = 1
       genomecov_2x(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_2x()")
  print("hg38_selected_genome =", hg38_selected_genome)
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("num_samples =", num_samples)
  print("num_test =", num_test)
  print("".join(["="] * 80)) # section-separating line

  # run bedtools genomecov bg on num_samples in parallel
  os.system("mkdir " + output_dir)
  process_multi_samples(genomecov_2x_per_bam, bam_list_file, output_dir, num_samples, num_test)

  # use {chr_name: array} to represent genome positions 
  genome_array = {}
  chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
  with open(hg38_selected_genome, 'r') as file_handle:
    for line in file_handle:
      chr_name, chr_length = line.split()
      chr_length = int(chr_length)
      genome_array[chr_name]  = np.zeros(chr_length, dtype='int16')
  ref_length = sum([len(x) for x in genome_array.values()])

  # summarize #samples with depth 2x at each genome position
  # limit the analysis to only chr1
  chr_name = "chr1"
  bed_file_list = [output_dir + x for x in os.listdir(output_dir) if ".genomecov_2x" in x]
  print("len(bed_file_list) =", len(bed_file_list))
  for bed_file in bed_file_list:
    with open(bed_file, 'r') as file_handle:
      for line in file_handle:
        if chr_name in line:
          chr_name, start, end, depth = line.split()
          start, end, depth = [int(x) for x in [start, end, depth]]
          genome_array[chr_name][start:end] += 1
        else:
          break
  np_file = output_dir + chr_name + "_depth_2x_distr.npy"
  np.save(np_file, genome_array[chr_name])
  #genome_array[chr_name] = np.load(np_file)
  print(np_file)

  # print histogram and summary
  print("chr_name =", chr_name)
  print("len(genome_array[chr_name]) =", len(genome_array[chr_name]))
  print("np.sum(genome_array[chr_name] >= 1) =", np.sum(genome_array[chr_name] >= 1))
  print("np.max(genome_array[chr_name]) =", np.max(genome_array[chr_name]))
  print("np.median(genome_array[chr_name]) =", np.median(genome_array[chr_name]))
  png_file = output_dir + chr_name + "_depth_2x_distr.png"
  fig, ax = pyplot.subplots()
  pyplot.hist(genome_array[chr_name], bins=range(1, 3*int(np.median(genome_array[chr_name]))),
              facecolor='mediumblue', alpha=1.0, rwidth=1.0)
  ax.set_ylabel("Number of genome positions")
  ax.set_xlabel("Number of bam files with 2 reads at a genome position")
  pyplot.savefig(png_file)
  print(png_file)

  # print bed file
  bed_file = output_dir + chr_name + "_depth_2x_distr.bed"
  nonzero_index = np.flatnonzero(genome_array[chr_name])
  current_start = nonzero_index[0]
  current_end = current_start
  current_depth = genome_array[chr_name][current_start]
  with open(bed_file, 'w') as file_handle:
    for index in nonzero_index[1:]:
      depth = genome_array[chr_name][index]
      if (current_end == index - 1 and current_depth == depth):
        current_end = index
      else:
        line = '\t'.join([chr_name, str(current_start), str(current_end+1), str(current_depth)])
        line += '\n'
        file_handle.write(line)
        current_start = index
        current_end = current_start
        current_depth = genome_array[chr_name][current_start]
    line = '\t'.join([chr_name, str(current_start), str(current_end+1), str(current_depth)])
    line += '\n'
    file_handle.write(line)
  print(bed_file)


def genomecov_bg_per_bam(io_tuple):
  """Run bedtools genomecov bg on 1 bam file.
     bedtools can only use 1 thread.
  """

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]

  bed_file = output_dir + bam_file_name + ".bg"
  command = ["bedtools genomecov"]
  command += ["-bg"]
  command += ["-ibam", bam_file]
  command += [">", bed_file]
  command = " ".join(command)
  os.system(command)


def genomecov_bg_depth_per_chr(io_tuple):
  """
  """

  chr_name, chr_length, output_dir = io_tuple
  chr_array = np.zeros(chr_length, dtype='int32')
  chr_name_len = len(chr_name)

  bed_file_list = [output_dir + x for x in os.listdir(output_dir) if ".bam.bg" in x]
  depth_1x = np.zeros(len(bed_file_list))
  depth_2plus = np.zeros(len(bed_file_list))
  for index, bed_file in enumerate(bed_file_list):
    with open(bed_file, 'r') as file_handle:
      for line in file_handle:
        name, start, end, depth = line.split()
        if chr_name == name:
          start, end, depth = [int(x) for x in [start, end, depth]]
          chr_array[start:end] += depth
          num_positions = end - start
          if depth == 1:
            depth_1x[index] += num_positions
          elif depth >= 2:
            depth_2plus[index] += num_positions
  depth_nonzero = depth_1x + depth_2plus
  
  return chr_name, chr_array, depth_1x, depth_2plus, depth_nonzero


def genomecov_bg_write_per_chr(io_tuple):
  """
  """

  chr_name, chr_bw_file, chr_array, hg38_selected_genome = io_tuple
  
  temp_bg_file = chr_bw_file + ".bg"
  with open(temp_bg_file, 'w') as file_handle:
    nonzero_index = np.flatnonzero(chr_array)
    current_start = nonzero_index[0]
    current_end = current_start
    current_depth = chr_array[current_start]
    for index in nonzero_index[1:]:
      depth = chr_array[index]
      if (index == current_end + 1 and depth == current_depth):
        current_end = index
      else:
        line = '\t'.join([chr_name, str(current_start), str(current_end+1), str(current_depth)])
        line += '\n'
        file_handle.write(line)
        current_start = index
        current_end = current_start
        current_depth = chr_array[current_start]
    line = '\t'.join([chr_name, str(current_start), str(current_end+1), str(current_depth)])
    line += '\n'
    file_handle.write(line)

  os.system("bedGraphToBigWig {0:s} {1:s} {2:s}".format(temp_bg_file, hg38_selected_genome, chr_bw_file))
  os.system("rm {0:s}".format(temp_bg_file))


def genomecov_bg(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test):
  """Run bedtools genomecov bg on a list of bam files.
     Calculate merged bed graph and sequencing depth.

     Usage:
       hg38_selected_genome = "/data/nh2tran/GeneSolutions/references/hg38_selected.genome"
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/genomecov_bg/"
       num_samples = 1
       num_test = 1
       genomecov_bg(hg38_selected_genome, bam_list_file, output_dir, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("genomecov_bg()")
  print("hg38_selected_genome =", hg38_selected_genome)
  print("bam_list_file =", bam_list_file)
  print("output_dir =", output_dir)
  print("num_samples =", num_samples)
  print("num_test =", num_test)
  print("".join(["="] * 80)) # section-separating line

  start_time = time.time()
  # run bedtools genomecov bg on num_samples in parallel
  os.system("mkdir " + output_dir)
  process_multi_samples(genomecov_bg_per_bam, bam_list_file, output_dir, num_samples, num_test)
  print("time.time() - start_time =", time.time() - start_time)
  print("".join(["="] * 80)) # section-separating line
  start_time = time.time()
  
  # aggregate depth along bed- and genome- dimentions
  # parallel over chr
  chr_list = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
  chr_length = []
  with open(hg38_selected_genome, 'r') as file_handle:
    for line in file_handle:
      _, length = line.split()
      chr_length.append(int(length))
  output_dir_list = [output_dir] * len(chr_list)
  argument_list = zip(chr_list, chr_length, output_dir_list)
  pool = multiprocessing.Pool(processes=num_samples)
  result_batch = pool.map(genomecov_bg_depth_per_chr, argument_list)
  pool.close()
  pool.join()
  genome_array = {}
  depth_1x, depth_2plus, depth_nonzero = 0, 0, 0
  for result in result_batch: # chr_name, chr_array, depth_1x, depth_2plus, depth_nonzero
    genome_array[result[0]] = result[1]
    depth_1x += result[2]
    depth_2plus += result[3]
    depth_nonzero += result[4]
  # make sure no overlapping reads left per bam
  ref_length = sum(chr_length)
  pct_nonzero = np.mean(depth_nonzero / ref_length)
  pct_1x = np.mean(depth_1x / ref_length)
  pct_2plus = np.mean(depth_2plus / ref_length)
  print("pct_nonzero per bam =", pct_nonzero)
  print("pct_1x =", pct_1x)
  print("pct_2plus =", pct_2plus)
  print("time.time() - start_time =", time.time() - start_time)
  print("".join(["="] * 80)) # section-separating line
  start_time = time.time()

  # calculate distribution of sequencing depth
  genomecov_hist_file = output_dir + "merged.genomecov_hist.full"
  with open(genomecov_hist_file, 'w') as file_handle:
    genome_hist = {}
    for chr_name in chr_list:
      depth, frequency = np.unique(genome_array[chr_name], return_counts=True)
      for x, y in zip(depth, frequency):
        line = '\t'.join([chr_name, str(x), str(y)])
        line += '\n'
        file_handle.write(line)
        if x in genome_hist:
          genome_hist[x] += y
        else:
          genome_hist[x] = y
    genome_depth = sorted(genome_hist.keys())
    for depth in genome_depth:
      line = '\t'.join(["genome", str(depth), str(genome_hist[depth])])
      line += '\n'
      file_handle.write(line)
  print("genomecov_hist_file =", genomecov_hist_file)
  print("time.time() - start_time =", time.time() - start_time)
  print("".join(["="] * 80)) # section-separating line
  start_time = time.time()

  # print bigwig file
  argument_list = [[chr_name, output_dir + "merged_" + chr_name + ".bw", genome_array[chr_name], hg38_selected_genome]
                   for chr_name in chr_list]
  pool = multiprocessing.Pool(processes=num_samples)
  pool.map(genomecov_bg_write_per_chr, argument_list)
  pool.close()
  pool.join()
  print("write merged_chr*.bw")
  print("time.time() - start_time =", time.time() - start_time)
  print("".join(["="] * 80)) # section-separating line
  start_time = time.time()


def filter_q30_1read(io_tuple):
  """Filter alignments with mapq >= 30.
     Sample 1 overlapping read (also remove secondary alignments).

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/bam.list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/filter_q30_1read/"
       num_samples = 1 # 1 threads
       num_test = 1
       filter_q30_1read(bam_list_file, output_dir, num_samples, num_test)
  """

  start_time = time.time()

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]
  prefix = bam_file_name[:-4]

  num_threads = 1

  # bam to sam, filter alignments with mapq >= 30
  temp_q30_sam = output_dir + prefix + ".q30.sam"
  command = ["samtools view -h"]
  command += ["-q", "30"]
  command += ["-F", "0x100"] # skip secondary alignments
  command += ["-@", str(num_threads)]
  command += [bam_file]
  command += ["-o", temp_q30_sam]
  command = " ".join(command)
  os.system(command)

  # use {chr_name: array} to represent genome regions, 0-based
  genome_array = {}
  with open("/data/nh2tran/GeneSolutions/references/hg38_selected.genome", 'r') as file_handle:
    for line in file_handle:
      chr_name, chr_length = line.split()
      chr_length = int(chr_length)
      genome_array[chr_name]  = np.zeros(chr_length, dtype=bool)
  # if a read falls into an empty region, write it and mark the region as True;
  # otherwise skip it
  temp_q30_1read_sam = output_dir + prefix + ".q30.1read.sam"
  with open(temp_q30_sam, 'r') as input_handle:
    with open(temp_q30_1read_sam, 'w') as output_handle:
      for line in input_handle:
        if line[0] == "@": # write header line
          output_handle.write(line)
        else:
          line_split = line.split()
          chr_name = line_split[2]
          if chr_name == "*": # unmapped read without coordinate, just write, don't mark anything
            output_handle.write(line)
          else:
            chr_pos = int(line_split[3]) - 1 # convert coordinate from 1-based to 0-based
            if np.any(genome_array[chr_name][chr_pos:chr_pos+75]): # skip overlapping read
              continue
            else: # empty region, write read and mark the region as True;
              output_handle.write(line)
              genome_array[chr_name][chr_pos:chr_pos+75] = True

  # sam to bam and index
  output_bam = output_dir + prefix + ".q30.1read.bam"
  command_view = ["samtools view"]
  command_view += ["-@", str(num_threads)]
  command_view += [temp_q30_1read_sam]
  command_view += ["-b -o", output_bam]
  command_view = " ".join(command_view)
  os.system(command_view)
  command_index = ["samtools index"]
  command_index += [output_bam]
  command_index = " ".join(command_index)
  os.system(command_index)

  # remove temp files
  os.system("rm " + temp_q30_sam)
  os.system("rm " + temp_q30_1read_sam)


def samtools_flagstat_per_bam(io_tuple):

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]

  num_threads = 4

  log_file = output_dir + bam_file_name + ".flagstat"
  command = ["samtools flagstat"]
  command += ["-@", str(num_threads)]
  command += [bam_file]
  command += [">", log_file]
  command = " ".join(command)
  os.system(command)


def samtools_flagstat(bam_list_file, output_dir, num_samples=1, num_test=1):
  """Summarize flagstat of a list bam files.

     Usage:
       bam_list_file = "/data/nh2tran/GeneSolutions/temp/step_3_filter_q20_1read.bam_list"
       output_dir = "/data/nh2tran/GeneSolutions/temp/step_5_flagstat/"
       num_samples = 1 # 4 threads
       num_test = 1
       samtools_flagstat(bam_list_file, output_dir, num_samples, num_test)
  """

  print("".join(["="] * 80)) # section-separating line
  print("samtools_flagstat()")
  print("bam_list_file =", bam_list_file)
  print("".join(["="] * 80)) # section-separating line

  os.system("mkdir " + output_dir)
  process_multi_samples(samtools_flagstat_per_bam, bam_list_file, output_dir, num_samples, num_test)

  log_file_list = [output_dir + x for x in os.listdir(output_dir) if ".bam.flagstat" in x]
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


def gatk_per_chromosome(args):

  bam_list_file, reference_file, chromosome = args

  command_list = ["gatk --java-options '-Xmx8G'"]

  # ~ command_list += ["HaplotypeCaller"]
  # Ploidy (number of chromosomes) per sample.
  # For pooled data, set to (Number of samples in each pool * Sample Ploidy).
  # Default value: 2. 
  # ~ command_list += ["--sample-ploidy", str(10)]

  command_list.append("Mutect2")
  command_list.append("--tumor-sample NIPT_700")

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
  #os.system("rm temp_chr*")


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


def draw_figure_2_venn():

  print("".join(["="] * 80)) # section-separating line
  print("draw_figure_2_venn()")
  print("".join(["="] * 80)) # section-separating line

  subsets = [597342,
             0,
             0,
             9947022,
             67153,
             5320214,
             7390020]
  nipt_snps =  8054515
  khv_snps =  12710234
  eas_snps =  22724409
  set_labels = ['A', 'B', 'C']
  venn_plot = venn3(subsets=subsets, set_labels=set_labels, set_colors=['r', 'b', 'g'])
  venn_plot.get_label_by_id('A').set_text('NIPT: {0:,} SNPs'.format(nipt_snps))
  venn_plot.get_label_by_id('B').set_text('KHV: {0:,} SNPs'.format(khv_snps))
  venn_plot.get_label_by_id('C').set_text('EAS: {0:,} SNPs'.format(eas_snps))
  venn_plot.get_label_by_id('A').set_color('r')
  venn_plot.get_label_by_id('A').set_alpha(1.0)
  venn_plot.get_label_by_id('B').set_color('b')
  venn_plot.get_label_by_id('B').set_alpha(1.0)
  venn_plot.get_label_by_id('C').set_color('g')
  venn_plot.get_label_by_id('C').set_alpha(1.0)
  venn_plot.get_label_by_id('100').set_text('{0:,}\n({1:.1%})'.format(subsets[0], subsets[0]/nipt_snps))
  venn_plot.get_label_by_id('101').set_text('{0:,}\n({1:.1%})'.format(subsets[4], subsets[4]/nipt_snps))
  venn_plot.get_label_by_id('111').set_text('{0:,}\n({1:.1%})'.format(subsets[6], subsets[6]/nipt_snps))
  venn_plot.get_label_by_id('011').set_text('')
  venn_plot.get_label_by_id('001').set_text('')
  pyplot.savefig("plot_figure_2_venn.png")
  pyplot.close()

  subsets = [664495,
             0,
             0,
             4599704,
             501004,
             720510,
             6889016]
  nipt_snps =  8054515
  khv_2pct_snps =  7609526
  khv_snps =  12710234
  set_labels = ['A', 'B', 'C']
  venn_plot = venn3(subsets=subsets, set_labels=set_labels, set_colors=['r', 'black', 'dodgerblue'])
  venn_plot.get_label_by_id('A').set_text('NIPT: {0:,} SNPs'.format(nipt_snps))
  venn_plot.get_label_by_id('B').set_text('KHV (AF >= 0.02):\n{0:,} SNPs'.format(khv_2pct_snps))
  venn_plot.get_label_by_id('C').set_text('KHV: {0:,} SNPs'.format(khv_snps))
  venn_plot.get_label_by_id('A').set_color('r')
  venn_plot.get_label_by_id('A').set_alpha(1.0)
  venn_plot.get_label_by_id('B').set_color('black')
  venn_plot.get_label_by_id('B').set_alpha(1.0)
  venn_plot.get_label_by_id('C').set_color('b')
  venn_plot.get_label_by_id('C').set_alpha(1.0)
  venn_plot.get_label_by_id('100').set_text('{0:,}\n({1:.1%})'.format(subsets[0], subsets[0]/nipt_snps))
  venn_plot.get_label_by_id('101').set_text('{0:,}\n({1:.1%})'.format(subsets[4], subsets[4]/nipt_snps))
  venn_plot.get_label_by_id('111').set_text('{0:,}\n({1:.1%})'.format(subsets[6], subsets[6]/nipt_snps))
  venn_plot.get_label_by_id('011').set_text('')
  venn_plot.get_label_by_id('001').set_text('')
  pyplot.savefig("plot_figure_2_venn_supp.png")
  pyplot.close()


def get_AF(isec_dir):
  """Get allele frequency from 1000 genomes or estimate from alt/ref read count.

     Usage:
       input_file =  "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/"
       sample_AF_file, population_AF_file = get_AF(isec_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("get_AF()")
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
      
      # read a variant site and its pair of (sample, population) AF
      num_sites = 0
      num_alleles = 0
      sample_AF_list = []
      population_AF_list = []
      while line_1 and line_2:
        num_sites += 1
        chr_1, position_1, _, ref_1, alt_1, _, _, info_1, _, sample_1 = line_1.split('\t')
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
        
        # sample allele depth (AD): ref_count,alt_count[,alt_count]
        sample_AD = sample_1.split(':')[1]
        sample_AD = sample_AD.split(',')
        sample_AD = [int(x) for x in sample_AD]
        sample_AD_total = float(sum(sample_AD))
        # exclude ref_count
        sample_AD = sample_AD[1:]
        # pair count to allele to avoid inconsistent order 
        sample_AD_dict = dict(zip(alt_1, sample_AD))
        # Mutect2 estimation of sample AF: alt_freq[,alt_freq]
        sample_AF = [x.split('=')[1] for x in info_1.split(';') if x[:8] == 'NIPT_AF=']
        sample_AF = sample_AF[0]
        sample_AF = sample_AF.split(',')
        sample_AF = [float(x) for x in sample_AF]
        sample_AF_dict = dict(zip(alt_1, sample_AF))
        
        # population AF: alt_freq[,alt_freq]
        population_AF = [x.split('=')[1] for x in info_2.split(';') if x[:3] == 'AF=']
        population_AF = population_AF[0]
        population_AF = population_AF.split(',')
        population_AF = [float(x) for x in population_AF]
        # pair freq to allele to avoid inconsistent order 
        population_AF_dict = dict(zip(alt_2, population_AF))
        
        # sample_AF is estimated by Mutect2 from alt/ref read count
        # population_AF is from 1000 genomes
        assert len(sample_AD) == len(population_AF)
        assert len(sample_AF) == len(population_AF)
          
        for allele in population_AF_dict:
          #sample_AF_list.append(sample_AD_dict[allele] / sample_AD_total)
          sample_AF_list.append(sample_AF_dict[allele])
          population_AF_list.append(population_AF_dict[allele])
          
        # read next variant pair
        line_1 = handle_1.readline()
        line_2 = handle_2.readline()

  print("num_sites =", num_sites)
  print("num_alleles =", num_alleles)
  print("len(sample_AF_list) =", len(sample_AF_list))
  print("len(population_AF_list) =", len(population_AF_list))
  sample_AF_file = isec_dir + "sample_AF_list.npy"
  population_AF_file = isec_dir + "population_AF_list.npy"
  np.save(sample_AF_file, sample_AF_list)
  np.save(population_AF_file, population_AF_list)
  return sample_AF_file, population_AF_file


def get_AF_EAS(isec_dir):
  """Get allele frequency from 1000 genomes or estimate from alt/ref read count.

     Usage:
       input_file =  "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/"
       sample_AF_file, population_AF_file = get_AF(isec_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("get_AF()")
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
      
      # read a variant site and its pair of (sample, population) AF
      num_sites = 0
      num_alleles = 0
      sample_AF_list = []
      population_AF_list = []
      EAS_AF_list = []
      while line_1 and line_2:
        num_sites += 1
        chr_1, position_1, _, ref_1, alt_1, _, _, info_1 = line_1.split('\t')
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
        
        # population AF: alt_freq[,alt_freq]
        info_1_split = info_1.split(';')
        sample_AF = [x.split('=')[1] for x in info_1_split if x[:3] == 'AF=']
        sample_AF = sample_AF[0]
        sample_AF = sample_AF.split(',')
        sample_AF = [float(x) for x in sample_AF]
        # pair freq to allele to avoid inconsistent order 
        sample_AF_dict = dict(zip(alt_2, sample_AF))

        # population AF: alt_freq[,alt_freq]
        info_2_split = info_2.split(';')
        population_AF = [x.split('=')[1] for x in info_2_split if x[:3] == 'AF=']
        population_AF = population_AF[0]
        population_AF = population_AF.split(',')
        population_AF = [float(x) for x in population_AF]
        # pair freq to allele to avoid inconsistent order 
        population_AF_dict = dict(zip(alt_2, population_AF))
        
        # sample_AF is estimated by Mutect2 from alt/ref read count
        # population_AF is from 1000 genomes
        assert len(sample_AF) == len(population_AF)
          
        for allele in population_AF_dict:
          sample_AF_list.append(sample_AF_dict[allele])
          population_AF_list.append(population_AF_dict[allele])
          
        # read next variant pair
        line_1 = handle_1.readline()
        line_2 = handle_2.readline()

  print("num_sites =", num_sites)
  print("num_alleles =", num_alleles)
  print("len(sample_AF_list) =", len(sample_AF_list))
  print("len(population_AF_list) =", len(population_AF_list))
  sample_AF_file = isec_dir + "sample_AF_list.npy"
  population_AF_file = isec_dir + "population_AF_list.npy"
  np.save(sample_AF_file, sample_AF_list)
  np.save(population_AF_file, population_AF_list)
  return sample_AF_file, population_AF_file


def compare_AF(sample_AF_file, population_AF_file, AF_max, scatter_plot_res):
  """Compare allele frequency from sample alt/ref read count and from 1000 genomes.

     Usage:
       sample_AF_file =   "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/sample_AF_list.npy"
       population_AF_file =   "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/population_AF_list.npy"
       AF_max = 1.0
       scatter_plot_res = 1
       compare_AF(sample_AF_file, population_AF_file, scatter_plot_res)
  """

  print("".join(["="] * 80)) # section-separating line
  print("compare_AF()")
  print("sample_AF_file =", sample_AF_file)
  print("population_AF_file =", population_AF_file)
  print("AF_max =", AF_max)
  print("scatter_plot_res =", scatter_plot_res)
  print("".join(["="] * 80)) # section-separating line

  sample_AF_list = np.load(sample_AF_file)
  population_AF_list = np.load(population_AF_file)
  threshold_index = population_AF_list <= AF_max
  sample_AF_list = sample_AF_list[threshold_index]
  population_AF_list = population_AF_list[threshold_index]
  print("sample_AF_list.shape =", sample_AF_list.shape)
  print("population_AF_list.shape =", population_AF_list.shape)
  corr = np.corrcoef(sample_AF_list, population_AF_list)[0, 1]
  error = sample_AF_list - population_AF_list
  ae = np.absolute(error)
  ape = ae / ((sample_AF_list + population_AF_list) / 2)
  mae = np.mean(ae)
  mape = np.mean(ape)
  median_ae = np.median(ae)
  median_ape = np.median(ape)
  pct_negative_error = np.sum(error < 0) / len(error)
  print("corr =", corr)
  print("mae =", mae)
  print("mape =", mape)
  print("median_ae =", median_ae)
  print("median_ape =", median_ape)
  print("pct_negative_error =", pct_negative_error)
  
  slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x=sample_AF_list, y=population_AF_list)
  print("slope =", slope)
  print("intercept =", intercept)
  print("r_value =", r_value)

  down_sampling = np.random.choice(population_AF_list.shape[0], size=population_AF_list.shape[0]//scatter_plot_res, replace=False)
  fig, ax = pyplot.subplots(figsize=(8, 8))
  pyplot.scatter(x=population_AF_list[down_sampling], y=sample_AF_list[down_sampling],
                 marker='.', s=0.1, color='salmon')
  ax.set_xlabel('KHV allele frequency', fontsize=18)
  ax.set_ylabel('NIPT allele frequency', fontsize=18)
  # ~ ax.set_xlim([0.0, AF_max])
  # ~ ax.set_ylim([0.0, AF_max])
  pyplot.text(0.0, 0.9, 'Pearson correlation = {0:.1%}'.format(corr), fontsize=18)
  ax.tick_params(labelsize=12)
  pyplot.savefig('plot_scatter_AF.png')
  pyplot.close()


def draw_overlap_venn():

  print("".join(["="] * 80)) # section-separating line
  print("draw_overlap_venn()")
  print("".join(["="] * 80)) # section-separating line

  records = [8632715, 7043477, 9452907]
  snps = [7610740, 6430459, 7655420]
  indels = [1043050, 614918, 1758977]
  khv, overlap, nipt = records
  khv_only = khv - overlap
  nipt_only = nipt - overlap
  venn_plot = venn2(subsets=(khv_only, nipt_only, overlap),
                    set_labels=['KHV (AF >= 2%)', 'NIPT_700_Mutect2'])
  venn_plot.get_label_by_id('10').set_text('{0:,}\n({1:.0%} KHV)'.format(khv_only, khv_only/khv))
  venn_plot.get_patch_by_id('10').set_color('red')
  venn_plot.get_patch_by_id('10').set_edgecolor('none')
  venn_plot.get_patch_by_id('10').set_alpha(0.5)
  venn_plot.get_label_by_id('11').set_text('{0:,}\n({1:.0%} KHV)'.format(overlap, overlap/khv))
  venn_plot.get_patch_by_id('11').set_color('blue')
  venn_plot.get_patch_by_id('11').set_edgecolor('none')
  venn_plot.get_patch_by_id('11').set_alpha(0.5)
  venn_plot.get_label_by_id('01').set_text('{0:,}'.format(nipt_only))
  pyplot.savefig("venn2.png")


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


def sample_2plus_reads(io_tuple):

  bam_file, output_dir = io_tuple
  bam_file_name = bam_file.split('/')[-1]
  prefix = bam_file_name[:-4]

  bed_file = output_dir + prefix + ".temp.depth_2plus.bed"
  command = ["bedtools genomecov"]
  command += ["-bg"]
  command += ["-ibam", bam_file]
  command += ["| grep -vw \"1$\""]
  command += [">", bed_file]
  command = " ".join(command)
  os.system(command)

  depth_1_bam = output_dir + prefix + ".temp.depth_1.bam"
  command = ["bedtools intersect"]
  command += ["-v"]
  command += ["-abam", bam_file]
  command += ["-b", bed_file]
  command += [">", depth_1_bam]
  command = " ".join(command)
  os.system(command)

  depth_2plus_bam = output_dir + prefix + ".temp.depth_2plus.bam"
  command = ["bedtools intersect"]
  command += ["-abam", bam_file]
  command += ["-b", bed_file]
  command += [">", depth_2plus_bam]
  command = " ".join(command)
  os.system(command)

  depth_2plus_downsampled_bam = output_dir + prefix + ".temp.depth_2plus_downsampled.bam"
  command = ["gatk DownsampleSam"]
  command += ["-I", depth_2plus_bam]
  command += ["-O", depth_2plus_downsampled_bam]
  command += ["-P", str(0.5)]
  command = " ".join(command)
  os.system(command)

  output_bam = output_dir + prefix + ".sampled_2plus.bam"
  command = ["samtools merge"]
  command += [output_bam]
  command += [depth_1_bam]
  command += [depth_2plus_downsampled_bam]
  command = " ".join(command)
  os.system(command)

  command = ["gatk BuildBamIndex"]
  command += ["-I", output_bam]
  command = " ".join(command)
  os.system(command)

  temp_files = output_dir + prefix + ".temp.*"
  command = ["rm", temp_files]
  command = " ".join(command)
  os.system(command)


def create_chr_mask(bed_file, hist_file, output_dir):
  """Create a mask of (mean +/- 3x sd) for each chromosome.

     Usage:
       bed_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.bedgraph"
       hist_file = "/data/nh2tran/GeneSolutions/temp/NIPT_700.q_10.bam.hist"
       output_dir = "/data/nh2tran/GeneSolutions/temp/create_chr_mask"
       create_chr_mask(bed_file, hist_file, output_dir)
  """

  print("".join(["="] * 80)) # section-separating line
  print("create_chr_mask()")
  print("bed_file =", bed_file)
  print("hist_file =", hist_file)
  print("output_dir =", output_dir)
  print("".join(["="] * 80)) # section-separating line
  
  assert not os.path.exists(output_dir), "output_dir exists " + output_dir
  os.makedirs(output_dir)

  chr_list = ["chr" + str(x) for x in range(1, 22+1)] + ["chrX", "chrY"]
  lower_bound = {}
  upper_bound = {}
  for chr_name in chr_list:
    median_depth, median_sd = genomecov_merge(hist_file, chr_name, False)
    lower_bound[chr_name] = median_depth - 3 * median_sd
    upper_bound[chr_name] = median_depth + 3 * median_sd
  chr_current = ""
  chr_file = ""
  chr_handle = None
  with open(bed_file, 'r') as file_handle:
    for line in file_handle:
      chr_name, start, end, depth = line.split()
      # if new chromosome, close the current bed file and open a new one
      if chr_name != chr_current:
        if chr_current != "":
          chr_handle.close()
          command = ["bedtools merge"]
          command += ["-i", chr_file]
          command += [">", chr_file + ".merged"]
          command = " ".join(command)
          os.system(command)
          command = "rm " + chr_file
          os.system(command)
        chr_file = output_dir + chr_name + ".bed"
        chr_handle = open(chr_file, 'w')
        chr_current = chr_name
        print(chr_name, lower_bound[chr_name], upper_bound[chr_name])
      depth = float(depth)
      if (depth >= lower_bound[chr_name] and depth <= upper_bound[chr_name]):
        chr_handle.write("{0}\t{1}\t{2}\n".format(chr_name, start, end))
    chr_handle.close()
    command = ["bedtools merge"]
    command += ["-i", chr_file]
    command += [">", chr_file + ".merged"]
    command = " ".join(command)
    os.system(command)
    command = "rm " + chr_file
    os.system(command)


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


def temp_ipynb():
  """Back-up ipython notebook
  """

  #from ngs_modules import *
  #%matplotlib inline

  # calculate #bam at depth-2plus genome positions 
  bed_file_list = [output_dir + x for x in os.listdir(output_dir) if "2plus.bed" in x]
  print("len(bed_file_list) =", len(bed_file_list))
  # use {chr_name: array} to represent genome positions 
  genome_array = {}
  chr_list = ['chr' + str(x) for x in range(1, 23)]
  chr_list += ['chrX', 'chrY']
  with open("hg38_selected.genome", 'r') as file_handle:
    for line in file_handle:
      chr_name, chr_length = line.split()
      chr_length = int(chr_length)
      genome_array[chr_name]  = np.zeros(chr_length, dtype='int16')
  ref_length = sum([len(x) for x in genome_array.values()])
  print("ref_length =", ref_length)
  for bed_file in bed_file_list:
    print(bed_file)
    with open(bed_file, 'r') as file_handle:
      for line in file_handle:
        chr_name, start, end, depth = line.split()
        start, end, depth = [int(x) for x in [start, end, depth]]
        if depth == 2:
          genome_array[chr_name][start:end] += 1
  os.system("mkdir genome_array")
  for chr_name in genome_array:
    np.save("genome_array/" + chr_name, genome_array[chr_name])

  genome_array_dir = "/data/nh2tran/GeneSolutions/temp/genome_array/"
  genome_array = {}
  chr_list = ['chr' + str(x) for x in range(1, 23)]
  chr_list += ['chrX', 'chrY']
  for chr_name in chr_list[:1]:
      genome_array[chr_name] = np.load(genome_array_dir + chr_name + ".npy")
      print(chr_name)
      print(len(genome_array[chr_name]))
      print(np.sum(genome_array[chr_name] >= 1))
      print(np.max(genome_array[chr_name]))
      print(np.median(genome_array[chr_name]))
      fig, ax = pyplot.subplots()
      pyplot.hist(genome_array[chr_name], bins=range(1, 51))
      ax.set_ylabel("Number of genome positions")
      ax.set_xlabel("Number of bam files with 2 reads at a genome position")
      pyplot.savefig("chr1_depth_2x_samples.png")
      break
      with open("chr1_depth_2x_samples.bed", 'w') as file_handle:
          for index in np.flatnonzero(genome_array[chr_name]):
              line = '\t'.join(['chr1', str(index), str(index+1), str(genome_array[chr_name][index])])
              line += '\n'
              file_handle.write(line)

  isec_dir = "/data/nh2tran/GeneSolutions/temp/isec.NIPT_700.Mutect2.khv.AF_0.02/"
  sample_AF_file, population_AF_file = get_AF(isec_dir)
  sample_AF_file = isec_dir + "sample_AF_list.npy"
  population_AF_file = isec_dir + "population_AF_list.npy"
  scatter_plot_res = 100
  compare_AF(sample_AF_file, population_AF_file, scatter_plot_res)


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

