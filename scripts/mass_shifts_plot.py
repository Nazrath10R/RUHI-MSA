#This creates a plot with the mass_shift that have been analyzed trought the RUHI-MSA

import sys
import pandas as pd
from ggplot import *

file_name = sys.argv[1]

df = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Duodenum_output.txt", sep="\t")

#print (df.head())

p = ggplot(aes(x='Peptide'), data=df) +\
    geom_histogram(binwidth=1) +\
    labs(y= "Frequency (log)", x = "Mass_shift") +\
    ggtitle("Duodenum mass shift frequency") +\
    scale_y_log()
#print (p)

p.save("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Duodenum_Mass_shifts.txt")




df1 = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Small_intestine_output.txt", sep="\t")

#print (df.head())

p1 = ggplot(aes(x='Peptide'), data=df1) +\
    geom_histogram(binwidth=1) +\
    labs(y= "Frequency (log)", x = "Mass_shift") +\
    ggtitle("Small intestine mass shift frequency") +\
    scale_y_log()
#print (p)

p1.save("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Small_intestine_Mass_shifts.txt")




df2 = pd.read_csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Spleen_output.txt", sep="\t")

#print (df.head())

p2 = ggplot(aes(x='Peptide'), data=df2) +\
    geom_histogram(binwidth=1) +\
    labs(y= "Frequency (log)", x = "Mass_shift") +\
    ggtitle(" Spleen mass shift frequency") +\
    scale_y_log()
#print (p)

p2.save("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Small_intestine_Mass_shifts.txt")




