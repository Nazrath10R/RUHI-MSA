#This creates a plot with the mass_shift that have been analyzed trought the RUHI-MSA

import sys
import pandas as pd
from ggplot import *

file_name = sys.argv[1]

df = pd.read_csv("./Results/" + file_name + "_output.txt", sep="\t")

#print (df.head())

p = ggplot(aes(x='Peptide'), data=df) +\
    geom_histogram(binwidth=1) +\
    labs(y= "Frequency (log)", x = "Mass_shift") +\
    ggtitle("Mass shift frequency") +\
    scale_y_log()
#print (p)

p.save("./Results/Myplots/"+file_name+"_Mass_Shifts")
