 # This function takes a psm.tsv file from a MSfragger "open" search and creates
 # the input file for RUHI-MSA

import sys
import pandas as pd

file_name = sys.argv[1]
output_name = file_name

try:
    df_MSFragger_output = pd.read_csv("./data/analysis_results/" + file_name + ".tsv", sep="\t")
except:
    print ("Check the location and/or name of your file. The file must be at RUHI-MSA/data/analysis_results/")
    exit(1)

#select specific columns
df_MSFragger_output = df_MSFragger_output[['Peptide', 'Calibrated Observed Mass',
       'Delta Mass', 'PeptideProphet Probability']]

#select only row where PeptideProphet is confident the peptide is correctly identified
df_MSFragger_output =  df_MSFragger_output[(df_MSFragger_output['PeptideProphet Probability'] > 0.99)]

df_MSFragger_output =  df_MSFragger_output[(df_MSFragger_output['Delta Mass'] != 0)]

df_MSFragger_output = df_MSFragger_output.drop(columns=['PeptideProphet Probability'])

df_MSFragger_output.to_csv("./input/"+output_name+".txt", sep="\t", index=False, header=False)

#print (df_MSFragger_output.head())
