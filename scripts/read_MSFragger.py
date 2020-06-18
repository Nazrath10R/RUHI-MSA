import pandas as pd
import sys

#alowing the user to add name of file to analyze
if (len(sys.argv) < 2):
    print ("This program takes at least 1 argument which is the name of the file")
    quit()
elif (len(sys.argv) < 3):
    file_name = sys.argv[1]
    output_name = file_name
    print ("You can add the name of the outputfile as the second argument! If not, the output is the same as the file_name")
elif (len(sys.argv) < 4):
    file_name = sys.argv[1]
    output_name = sys.argv[2]
else:
    print ("This program at least takes 1 argument")
    quit()

try:
    df_MSFragger_output = pd.read_csv("../data/analysis_results/" + file_name, sep="\t")
except:
    print ("Check the location and/or name of your file. The file path should be is ../data/analysis_results/")
    exit(1)

#select specific columns
df_MSFragger_output = df_MSFragger_output[['Peptide', 'Observed M/Z',
       'Adjusted Delta Mass', 'PeptideProphet Probability', 'Assigned Modifications']]

#select only row where PeptideProphet is confident the peptide is correctly identified
df_MSFragger_output =  df_MSFragger_output[(df_MSFragger_output['PeptideProphet Probability'] > 0.99)]
#remove peptides where Modficications have been assigned

df_MSFragger_output =  df_MSFragger_output[(df_MSFragger_output['Adjusted Delta Mass'] != 0)]

df_MSFragger_output =  df_MSFragger_output[(df_MSFragger_output['Assigned Modifications'].isnull())]

df_MSFragger_output = df_MSFragger_output.drop(columns=['PeptideProphet Probability', 'Assigned Modifications'])

#df_MSFragger_output = df_MSFragger_output.drop_duplicates(keep='first', inplace=False)

df_MSFragger_output.to_csv("../input/"+output_name+".txt", sep="\t", index=False, header=False)

#print (df_MSFragger_output.head())
