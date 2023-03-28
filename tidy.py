#!/usr/bin/env python
import pandas as pd
import re

def parse_csv(excel_file):
  df = pd.read_excel(excel_file, header=0)

  # Remove uneccessary columns from bigg data
  if 'bigg' in excel_file:
    df = df.drop(labels=['database_links', 'model_list', 'old_bigg_ids'], axis=1)

  if excel_file == 'reactions.XLSX':
    # Fix formating of GPR column
    df['Gene to Protein to reaction association'] = df['Gene to Protein to reaction association'].map(reformat_text)
    df.insert(1, 'bigg_id', [None] * df.shape[0])
    df.insert(2, 'bigg_equation', [None] * df.shape[0])

  return df


def reformat_text(gpr_string):
  parens_list = []

  if gpr_string == '-':
    return '-'
  else:
    # Grab all text in ()
    reg_exp_parens = '\(.*?\)'
    parens_list = re.findall(reg_exp_parens, gpr_string)
    
    # Get rid of whitespace after ( and before )
    for i,paren in enumerate(parens_list):
      paren = re.sub('\s*\)', ')', paren)
      paren = re.sub('\(\s*', '(', paren)
      parens_list[i] = paren

    if len(parens_list) == 1:
      return parens_list[0]
    else:
      return ' or '.join(parens_list)

def format_exchange_items(exchange_string):
  string_list = exchange_string.split(' ')
  if string_list[0] == "exchange":
    new_string = string_list[-1] + " " + string_list[0]
    return new_string.lower()
  else:
    return None


def add_as_much_bigg_as_possible(react_df, bigg_df):
  count = 0
  for index,react_row in react_df.iterrows():
    count += 1
    if count % 5 == 0:
      left = react_df.shape[0] - count
      percent = left / react_df.shape[0]
      print('5 more completed! {}% complete.'.format(percent))

    exchange_name = format_exchange_items(react_row['Name'])
    
    for index, bigg_row in bigg_df.iterrows():
      if exchange_name != None and str(bigg_row['name']).lower() == exchange_name:
        react_row['bigg_id'] = bigg_row['bigg_id']
        react_row['bigg_equation'] = bigg_row['reaction_string']
        break
      elif react_row['Name'].lower() == str(bigg_row['name']).lower():
        react_row['bigg_id'] = bigg_row['bigg_id']
        react_row['bigg_equation'] = bigg_row['reaction_string']
        break
      else:
        pass

def main():
  metabolites_file = 'metabolites.XLSX'
  bigg_metabolites_file = 'bigg_models_metabolites.xlsx'
  reactions_file = 'reactions.XLSX'
  bigg_reactions_file = 'bigg_models_reactions.xlsx'

  #bigg_met_df = parse_csv(bigg_metabolites_file)
  bigg_reac_df = parse_csv(bigg_reactions_file)
  #met_df = parse_csv(metabolites_file)
  reac_df = parse_csv(reactions_file)
  add_as_much_bigg_as_possible(reac_df, bigg_reac_df)
  print(reac_df)
  reac_df.to_csv('test_auto_reactions.csv')


if __name__ == "__main__":
  main()