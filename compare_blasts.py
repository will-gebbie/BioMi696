#!/usr/bin/env python
import subprocess
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

sm22 = "/home/will_gebbie/blast_compare/sm22_prokka.faa"
xml_files = []

def file_to_string(filename):
  f = open(filename, "r")
  text = f.read()
  f.close()
  return text


def blast_align_cmd(s1, s2):
  # Gets file name without extension
  s1_name = s1.split('/')[-1].split('.')[0]
  s2_name = s2.split('/')[-1].split('.')[0]

  out_name = "{}_compare_with_{}.xml".format(s1_name, s2_name)

  return NcbiblastpCommandline(query=s1, subject=s2, out=out_name, outfmt=5, evalue=0.0001, qcov_hsp_perc=85, num_threads=8)


def convert_to_str_list(cmd):
  subprocess_input = cmd.split()
  return subprocess_input


def blast_all(prot_list):
  count = 0
  for prot in prot_list:
    cmd = convert_to_str_list(str(blast_align_cmd(sm22, prot)))
    count += 1
    print("Started blastp #{}!".format(count))
    subprocess.run(cmd)
    print("Finished blastp #{}!".format(count))
  
  print("All blastp's have finished!")


def parse_xml(xml_file):
  
  f = open(xml_file, 'r')
  write_filename = xml_file.split('_compare_with_')
  write_filename = "{}_genes_mapped_to_{}_genes.csv".format(write_filename[0], write_filename[1].split('.')[0])
  write_file = open(write_filename, 'w')

  add_header_to_file(write_file)

  count = 0
  xml_iter = NCBIXML.parse(f)
  
  for item in xml_iter:
    count += 1
    align_info = get_alignment_info(item)
    query_info = get_query_info(item)
    add_query_align_to_file(write_file, count, query_info, align_info)
  
  f.close()
  write_file.close()


def get_alignment_info(blast_record):
  record_info = {"Hit_num": [], "Hit_def": [], "Hit_len": [], "Hsp_score": [], "Hsp_evalue": [], "Hsp_align-len": []}
  c = 0
  for info in blast_record.alignments:
    
    c += 1
    
    record_info["Hit_num"].append(c)
    record_info["Hit_def"].append(info.hit_def)
    record_info["Hit_len"].append(info.length)
    
    for hsp_info in info.hsps:
      record_info["Hsp_score"].append(hsp_info.score)
      record_info["Hsp_evalue"].append(hsp_info.expect)
      record_info["Hsp_align-len"].append(hsp_info.align_length)
  
  # No Hits for this record?
  if c == 0:
    return {}
  else:
    return record_info

def get_query_info(blast_record):
  return blast_record.query


def add_query_align_to_file(file_handle, iter_num, query_data, align_info):
  if bool(align_info) is True:
    hit_nums = align_info["Hit_num"]
    hit_defs = align_info["Hit_def"]
    hit_lens = align_info["Hit_len"]
    hsp_scores = align_info["Hsp_score"]
    hsp_evalues = align_info["Hsp_evalue"]
    hsp_align_lens = align_info["Hsp_align-len"]
    num_hits = len(hit_nums)
    for i in range(num_hits):
      file_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(iter_num, query_data, hit_nums[i], hit_defs[i], hit_lens[i], hsp_scores[i], hsp_evalues[i], hsp_align_lens[i]))


def add_header_to_file(file_handle):
  header_string = "Iteration_iter-num\tIteration_query-def\tHit_num\tHit_def\tHit_len\tHsp_score\tHsp_evalue\tHsp_align-len\n"
  file_handle.write(header_string)
  

def main():
  
  methylo = "/home/will_gebbie/blast_compare/Methylobacterium.faa"
  bath = "/home/will_gebbie/blast_compare/McBath_work.fasta"
  nostoc = "/home/will_gebbie/blast_compare/Nostoc7120_ready.fasta"
  clostridium = "/home/will_gebbie/blast_compare/clostridium.faa"
  tricho = "/home/will_gebbie/blast_compare/trichosporium_ob3b.faa"
  rhodo = "/home/will_gebbie/blast_compare/rhodo.faa"
  extorquens_am1 = "/home/will_gebbie/blast_compare/extorquens_am1.fasta"

  compare_list = [methylo, bath, nostoc, clostridium, tricho, rhodo, extorquens_am1]

  # Uncomment if need to blast
  # blast_all(compare_list)

  dirname = "/home/will_gebbie/blast_compare/"
  ext = 'xml'

  for file in os.listdir(dirname):
    if file.endswith(ext):
      print("Collecting important data from {}...".format(file))
      parse_xml(file)
      print("Done")


main()

