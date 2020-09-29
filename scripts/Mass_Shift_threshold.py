#This creates a plot with the mass_shift that have been analyzed
# by the RUHI-MSA excluding the small Mass_Shifts (-1.5 to 1.5)

import plotly.express as px
#import plotly
import pandas as pd
import os
##---------------File 1 -----------------------------------------------------------------##

df = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Spleen_output.txt",sep="\t", index_col=False)
#df = df.replace("Unknown:177", "NaN")
#print (df.head(n=8))
df = df.iloc[:,[5,6,7,8,9]]
df = df[~df.iloc[:, 0].between(-1.5, 1.5, inclusive=False)]
#df.reset_index()
#print (df.columns)
#counting PTMs
df= df.apply(pd.value_counts)

df = df.fillna(0)
df = df[df.index != "Unknown:177"]
print (df)
df['Sum'] = df.sum(axis=1)
df= df.drop(columns= ['PTM1', 'PTM2', 'PTM3', 'PTM4', 'PTM5'])

df["PTMs"] = df.index
df["Sample"] = "Spleen"
#print (df.head())


##---------------File 2 -----------------------------------------------------------------##
df1 = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Duodenum_output.txt",sep="\t", index_col=False)


df1 = df1.iloc[:,[5,6,7,8,9]]
#df1.reset_index()
#print (df1.columns)
#counting PTMs
df1= df1.apply(pd.value_counts)

df1 = df1.fillna(0)