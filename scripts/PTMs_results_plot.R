library(ggplot2)

df <- read.csv("JOHNPP_l1_output.txt", sep = "\t", row.names = NULL)

#define column names
colnames(df)<- c('Peptide', ' Mass_shift', 'Peptide.Mass', 'Mass.of.PTMs', 'Score', 'PTM.1', 'PTM.2','PTM.3','PTM4', 'PTM5')

df <- unique( df[ , 1:10 ] ) #removing duplicates

#define variable as sring
df$PTM.1 <- as.character(df$PTM.1)
df$PTM.2 <- as.character(df$PTM.2)
df$PTM.3 <- as.character(df$PTM.3)
df$PTM4 <- as.character(df$PTM4)
df$PTM5 <- as.character(df$PTM5)
#make dataframe with only PTM results
frequency_df <- as.data.frame(table(c(df$PTM.1, df$PTM.2, df$PTM.3, df$PTM4, df$PTM5)))


#organize in decreasing order
frequency_df <- frequency_df[order(frequency_df$Freq, decreasing = TRUE),]

#removing rows that have empty Var1
frequency_df <- frequency_df[!(frequency_df$Var1==""),]

#make it as factor to be able to show the plot in decreasing order
frequency_df$Var1 <- factor(frequency_df$Var1, levels = frequency_df$Var1[order(frequency_df$Freq, decreasing = TRUE)])

##need to remove PTM.1, 2,3 etc and firs row

#plot
ggplot(frequency_df, aes(x=Var1, y=Freq)) +
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs (x ="PTMs", y = "Frequency (log10)") + ggtitle("JOHNPP_l1_output.txt") + scale_y_log10()
