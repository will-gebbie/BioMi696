#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

EVAL_INDEX = 0
SCORE_INDEX = 1
ALIGN_LEN_INDEX = 2

def add_histogram(ax, data, datatype, title):
  hist = data
  ax.set(title='{}'.format(title), xlabel="Value", ylabel="Frequency")

  if datatype == 'e-value':
    ax.hist(hist, bins=20, range=(0, 1*(10**-60)), edgecolor='black', linewidth=1.2)
    ax.set_ylim(0, 3000)
  else:
    ax.hist(hist, bins=20, edgecolor='black', linewidth=1.2)
  
  if datatype == 'score':
    ax.set_ylim(0,1000)
  elif datatype == 'align-len':
    ax.set_ylim(0,2000)
  
  #fig.savefig('{}_vs_{}_{}_histogram.png'.format(names[0],names[1], datatype))


def get_data(dataframe):
  e_val_data = dataframe['Hsp_evalue'].values
  align_len_data = dataframe['Hsp_align-len'].values
  score_data = dataframe['Hsp_score'].values

  return [e_val_data, score_data, align_len_data]


def parse_csv(csv_file, hist_data_dict):
  # Text processing
  names = csv_file.split('/')
  names = names[-1].split('_genes_mapped_to_')
  name1 = names[0]
  name2 = names[-1].split('_genes.csv')[0]
  names = [name1, name2]


  df = pd.read_csv(csv_file, sep='\t', header=0)
  dict_key_name = "{} vs {}".format(names[0], names[1])
  hist_data_dict[dict_key_name] = get_data(df)


def create_subplots(data, dtype):
  datatype = list(data.keys())[0]
  values = data[datatype]
  titles = data['comparison_names']
  
  fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(20,15))

  for title, value, ax in zip(titles, values, axs.flatten()):
    add_histogram(ax, value, datatype, title)

  fig.delaxes(axs.flatten()[-1])
  fig.subplots_adjust(wspace=1, hspace=None)
  fig.suptitle('Comparison of {}s'.format(dtype), fontsize=12)
  fig.savefig('{} Comparisons.png'.format(dtype))
  plt.close(fig)


def populate_hist_lists(hist_data_dict, align_list, eval_list, score_list):
  for key,value in hist_data_dict.items():
    align_list['align-len'].append(value[ALIGN_LEN_INDEX])
    align_list['comparison_names'].append(key)

    eval_list['e-value'].append(value[EVAL_INDEX])
    eval_list['comparison_names'].append(key)

    score_list['score'].append(value[SCORE_INDEX]) 
    score_list['comparison_names'].append(key)


def main():

  dirname = "/home/will_gebbie/blast_compare/"
  ext = ('csv', 'png')

  hist_data_dict = {}

  for file in os.listdir(dirname):
    if file.endswith(ext[0]):
      print("Creating relevant histograms from {}...".format(file))
      parse_csv(file, hist_data_dict)
      print("Done")

  

  align_hist_data = {'align-len': [], 'comparison_names': []}
  e_val_data = {'e-value': [], 'comparison_names': []}
  score_data = {'score': [], 'comparison_names': []}

  populate_hist_lists(hist_data_dict, align_hist_data, e_val_data, score_data)
  
  create_subplots(align_hist_data, "Align-Length")
  create_subplots(e_val_data, "E Value")
  create_subplots(score_data, "Blast Score")
      
main()