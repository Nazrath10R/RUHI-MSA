import plotly.express as px
#import plotly
import pandas as pd
import os
##---------------File 1 -----------------------------------------------------------------##

df = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Spleen_output.txt",sep="\t", index_col=False)
df = df.iloc[:,[5,6,7,8,9]]
#df.reset_index()
#print (df.columns)
#counting PTMs
df= df.apply(pd.value_counts)

df = df.fillna(0)

df['Sum'] = df.sum(axis=1)
df= df.drop(columns= ['PTM1', 'PTM2', 'PTM3', 'PTM4', 'PTM5'])

df["PTMs"] = df.index
df["Sample"] = "Spleen"
#print (df.head())


##---------------File 2 -----------------------------------------------------------------##
df1 = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Duodenum_output.txt",sep="\t", index_col=False)
df1 = df1.replace("Unknown:177", 0)
df1 = df1.iloc[:,[5,6,7,8,9]]
#df1.reset_index()
#print (df1.columns)
#counting PTMs
df1= df1.apply(pd.value_counts)

df1 = df1.fillna(0)


df1['Sum'] = df1.sum(axis=1)
df1= df1.drop(columns= ['PTM1', 'PTM2', 'PTM3', 'PTM4', 'PTM5'])

df1["PTMs"] = df1.index
df1["Sample"] = "Duodenum"
#print (df1.head())

##---------------File 3 -----------------------------------------------------------------##
df2 = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Small_intestine_output.txt",sep="\t", index_col=False)
df2 = df2.iloc[:,[5,6,7,8,9]]
#df1.reset_index()
#print (df1.columns)
#counting PTMs
df2= df2.apply(pd.value_counts)

df2 = df2.fillna(0)

df2['Sum'] = df2.sum(axis=1)
df2= df2.drop(columns= ['PTM1', 'PTM2', 'PTM3', 'PTM4', 'PTM5'])

df2["PTMs"] = df2.index
df2["Sample"] = "Small intestine"
#print (df1.head())

##------------------Joining dataframes -----------------------------------------##
result = pd.concat([df,df1, df2])
#print (result.head())

##------------------Producing plot -----------------------------------------##


fig = px.bar(result, x="PTMs", y="Sum", color='Sample')
fig.update_layout(barmode='group', xaxis={'categoryorder':'total descending'}, yaxis_type="log")
#fig.show()
#plotly.offline.plot(fig, filename='./comparing_PTMs.html')
fig.write_html("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/comparing_3_tissues.html")
#fig.write_image("Users/pedrocardoso/Documents/bargraph.png")

print("FInished")





