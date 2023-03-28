#!/usr/bin/env python
import pandas as pd
import re

def parse_csv(excel_file):
  df = pd.read_excel(excel_file, header=0)

  # Remove uneccessary columns from bigg data
  if excel_file == 'bigg_models_metabolites.xlsx':
    df = df.drop(labels=['bigg_id', 'model_list', 'old_bigg_ids', 'database_links'], axis=1)

  elif excel_file == 'metabolites.XLSX' or excel_file == 'mapped_metabolites.xlsx':
    pass

  elif excel_file == 'mapped_bigg_reactions.xlsx':
    df = df[['bigg_equation', 'Equation']]

  return df


def match_mets(met_set, mapped_reactions, mets_df):
  f = open('persistent_met_set.txt', 'a')
  f2 = open('matched_mets.txt', 'a')
  for index, row in mapped_reactions.iterrows():
    curr_bigg_eq = str(row['bigg_equation'])
    curr_model_eq = str(row['Equation'])

    bigg_eq_mets, curr_model_mets = create_met_lists(curr_bigg_eq, curr_model_eq)
    
    for met in bigg_eq_mets:
      if met not in met_set:
        if met == 'nan':
          break
        met_set.add(met)
        model_met = manual_map(met, curr_model_mets)
        f.write(str(met) +'\n')
        f2.write(str(model_met) + '\t' + str(met) + '\n')
        
  f.close()
  f2.close()
    

def create_met_lists(eq1, eq2):
  comps = ['_e', '_c', '_p', '[e]', '[p]']
  mutipliers = ['0.5', '2.0', '2', '3.0', '3', '4.0', '4', '5.0', '5', '6.0', '6',
                '7.0', '7', '8.0', '8', '10.0','10']
  eq1list = []
  eq2list = []

  eq1 = eq1.replace('<->', '')
  eq1 = eq1.replace(' + ', ' ')
  eq2 = eq2.replace('=', '')
  eq2 = eq2.replace(' + ', ' ')

  # remove all compartment identifers
  for c in comps:
    if c in eq1:
      eq1 = eq1.replace(c, '')
    if c in eq2:
      eq2 = eq2.replace(c, '')

  # create list of metabolites and remove all whitespace
  eq1list = eq1.split()
  eq2list = eq2.split()

  # remove all numbers (multipliers of compounds in reactions)
  for m in mutipliers:
    if m in eq1list:
      eq1list = [x for x in eq1list if x != m]
    if m in eq2list:
      eq2list = [x for x in eq2list if x != m]

  return eq1list, eq2list


def manual_map(bigg_met, pos_model_mets):
  if bigg_met == 'nan' or bigg_met == 'nh4':
    return
  print("Which metabolite does this bigg_id match: {}".format(bigg_met))
  for i in range(len(pos_model_mets)):
    print(str(i) + ":", pos_model_mets[i], end='; ')

  print()
  index = -1
  while True:
    try:
      index = input("Enter # of correct metabolite: ")
      if int(index) in range(len(pos_model_mets)):
        break
      else:
        print("That is not a valid index, try again...")
    except ValueError:
      print("That is not a valid index, try again...")

  return pos_model_mets[int(index)]


def add_mets_to_file(match_mets, met_df):

  for index, row in match_mets.iterrows():
    model_met = str(row['Names'])
    bigg_id = str(row['Bigg'])
    print(model_met)
    # Find the index in the dataframe where the current metabolite is located
    if model_met != 'None':
      index_of_model_met = met_df['Name'].loc[lambda x: x==model_met].index[0]
      met_df.at[index_of_model_met, 'Universal Bigg Id'] = bigg_id


def main():
  metabolites_file = 'metabolites.XLSX'
  bigg_metabolites_file = 'bigg_models_metabolites.xlsx'
  mapped_reactions_file = 'mapped_bigg_reactions.xlsx'
  mapped_metabolites_file = 'mapped_metabolites.xlsx'

  # bigg_met_df = parse_csv(bigg_metabolites_file)
  # mapped_reac_df = parse_csv(mapped_reactions_file)
  met_df = parse_csv(metabolites_file)
  mapped_mets_df = parse_csv(mapped_metabolites_file)

  print()
  met_set = set()
  try:
    met_set = set(line.strip() for line in open('persistent_met_set.txt'))
  except FileNotFoundError:
    print('Creating persistent_met_set.txt...')

  # match_mets(met_set, mapped_reac_df, met_df)
  add_mets_to_file(mapped_mets_df, met_df)

  # Make new file
  met_df.to_csv('test_map_met.csv', index=False)
  

  


if __name__ == "__main__":
  main()


