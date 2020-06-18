import pandas as pd
from ggplot import *


df = pd.read_csv("/Users/pedrocardoso/Documents/RUHI-MSA/results/JOHNPP_l1_output.txt", sep="\t")

print (df.head())

p = ggplot(aes(x='Peptide'), data=df) +\
    geom_histogram(binwidth=0.01) +\
    labs(y= "Frequency (log)", x = "Mass shift") +\
    ggtitle("Mass shift frequency") +\
    scale_y_log()
#print (p)

p.save('/Users/pedrocardoso/Documents/RUHI-MSA/results/Myplots/JOHNPP_l1_mass_shifts.png')
