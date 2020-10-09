# This function concatenates files
 # the input file for RUHI-MSA

import sys
import pandas as pd

file_name = sys.argv[1]
output_name = file_name


file_name = "./Results/" + file_name




mydf = pd.read_csv(file_name + "_OUT.txt", sep="\t", index_col = False)
#print (mydf.head())
mydf = mydf.drop_duplicates()
#print (mydf.head())

mydf.to_csv(file_name + "_OUTPUT.txt",sep="\t", index = False)
# df1 = pandas.read_csv('path1')
# df2 = pandas.read_csv('path2', skiprows=1)
# df3 = pandas.read_csv('path3', skiprows=1)