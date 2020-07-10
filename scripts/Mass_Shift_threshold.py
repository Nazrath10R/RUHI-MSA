#This creates a plot with the mass_shift that have been analyzed
# by the RUHI-MSA excluding the small Mass_Shifts (-1.5 to 1.5)

import sys
import pandas as pd
from ggplot import *

file_name = sys.argv[1]

df = pd.read_csv("./Results/" + file_name + "_output.txt", sep="\t")
#print (df)
df = df[~df.iloc[:, 0].between(-1.5, 1.5, inclusive=False)]


p = ggplot(aes(x='Peptide'), data=df) +\
    geom_histogram(binwidth=0.5) +\
    labs(y= "Frequency (log)", x = "Mass_shift") +\
    ggtitle("Mass shift frequency") +\
    scale_y_log()
#print (p)

p.save("./Results/Myplots/"+file_name+"_Mass_Shifts_threshold")
