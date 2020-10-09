import plotly.express as px
#import plotly
import pandas as pd
import os
import sys

file_name = sys.argv[1]

df = pd.read_csv("./Results/" + file_name + "_OUTPUT.txt",sep="\t", index_col=False)
#print (df.head())
df = df.iloc[:,[5,6,7,8,9]]
#df.reset_index()
#print (df.columns)
#counting PTMs
df= df.apply(pd.value_counts)

df = df.fillna(0)

df['Sum'] = df.sum(axis=1)
df= df.drop(columns= ['PTM1', 'PTM2', 'PTM3', 'PTM4', 'PTM5'])

df["PTMs"] = df.index
df["Sample"] = file_name
#print (df.head())


##------------------Producing plot -----------------------------------------##


fig = px.bar(df, x="PTMs", y="Sum", color='Sample')
fig.update_layout(barmode='group', xaxis={'categoryorder':'total descending'}, yaxis_type="log")
#fig.show()
#plotly.offline.plot(fig, filename='./comparing_PTMs.html')
fig.write_html("./Results/Myplots/" + file_name + "_PTMs_PROFILE.html")
#fig.write_image("Users/pedrocardoso/Documents/bargraph.png")

#print("FInished")