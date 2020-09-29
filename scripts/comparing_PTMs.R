library(ggplot2)
library(dplyr)
install.packages("reshape2")
#args = commandArgs(trailingOnly=TRUE)

#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#} else if (length(args)==1) {
  
#}
#file <- paste("./Results/", args[1], "_output.txt", sep = "", collapse = NULL)
#file_output <- paste("./Results/Myplots/", args[1], "_PTMs.pdf", sep = "", collapse = NULL)


df_Small <- read.csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Small_intestine_output.txt",
                  sep = "\t", row.names = NULL)
df_Spleen <- read.csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Spleen_output.txt",
                sep = "\t", row.names = NULL)
df_Duodenum <- read.csv("/Users/pedrocardoso/Documents/PTMs_project/Sample_analysis/Tissue_results/Duodenum_output.txt",
                sep = "\t", row.names = NULL)

#define column names
colnames(df_Small)<- c('Peptide', ' Mass_shift', 'Peptide.Mass', 'Mass.of.PTMs', 'Score', 'PTM.1', 'PTM.2','PTM.3','PTM4', 'PTM5')
colnames(df_Spleen)<- c('Peptide', ' Mass_shift', 'Peptide.Mass', 'Mass.of.PTMs', 'Score', 'PTM.1', 'PTM.2','PTM.3','PTM4', 'PTM5')
colnames(df_Duodenum)<- c('Peptide', ' Mass_shift', 'Peptide.Mass', 'Mass.of.PTMs', 'Score', 'PTM.1', 'PTM.2','PTM.3','PTM4', 'PTM5')

#define variable as sring
#df$PTM.1 <- as.character(df$PTM.1)
#df$PTM.2 <- as.character(df$PTM.2)
#df$PTM.3 <- as.character(df$PTM.3)
#df$PTM4 <- as.character(df$PTM4)
#df$PTM5 <- as.character(df$PTM5)
#make dataframe with only PTM results
frequency_Small <- as.data.frame(table(c(df_Small$PTM.1, df_Small$PTM.2, df_Small$PTM.3, df_Small$PTM4, df_Small$PTM5)))

frequency_Spleen <- as.data.frame(table(c(df_Spleen$PTM.1, df_Spleen$PTM.2, df_Spleen$PTM.3, df_Spleen$PTM4, df_Spleen$PTM5)))

frequency_Duodenum <- as.data.frame(table(c(df_Duodenum$PTM.1, df_Duodenum$PTM.2, df_Duodenum$PTM.3, df_Duodenum$PTM4, df_Duodenum$PTM5)))


#organize in decreasing order
#frequency_df <- frequency_df[order(frequency_df$Freq, decreasing = TRUE),]
#frequency_df1 <- frequency_df1[order(frequency_df1$Freq, decreasing = TRUE),]

#removing rows that have empty Var1
frequency_Small <- frequency_Small[!(frequency_Small$Var1==""),]
frequency_Small<- frequency_Small[!(frequency_Small$Var1=="Unknown:177"),]
frequency_Small <- frequency_Small[!(frequency_Small$Var1=="#NAME!"),]
frequency_Spleen <- frequency_Spleen[!(frequency_Spleen$Var1==""),]
frequency_Spleen <- frequency_Spleen[!(frequency_Spleen$Var1=="Unknown:177"),]
frequency_Spleen <- frequency_Spleen[!(frequency_Spleen$Var1=="#NAME!"),]
frequency_Duodenum <- frequency_Duodenum[!(frequency_Duodenum$Var1==""),]
frequency_Duodenum <- frequency_Duodenum[!(frequency_Duodenum$Var1=="Unknown:177"),]
frequency_Duodenum <- frequency_Duodenum[!(frequency_Duodenum$Var1=="#NAME!"),]

#removing PTMs that show up less than 3 times
frequency_Small <- frequency_Small[!(frequency_Small$Freq < 1000),]
frequency_Spleen <- frequency_Spleen[!(frequency_Spleen$Freq < 1000),]
frequency_Duodenum <- frequency_Duodenum[!(frequency_Duodenum$Freq < 1000),]
#Create 1 column with name of sample
frequency_Small$Sample <- "Small intestine"
frequency_Spleen$Sample <- "Spleen"
frequency_Duodenum$Sample <- "Duodenum"

#joining the two dataframes
joined_df <- bind_rows(frequency_Small, frequency_Spleen, frequency_Duodenum)

#ordering
#joined_df <- joined_df[order(joined_df$Freq, decreasing = TRUE),]
#joined_df <- joined_df[!(joined_df$Freq < 5),] ##this is removing values from only one sample if one the sample has a freq > than 5

#make it as factor to be able to show the plot in decreasing order
#joined_df$Var1 <- factor(joined_df$Var1, levels = joined_df$Var1[order(joined_df$Freq, decreasing = TRUE)])

##need to remove PTM.1, 2,3 etc and firs row



#reshaping data format (for PCA)
new_df <- reshape(joined_df, idvar = "Var1", timevar = "Sample", direction = "wide")
#converting NA to 0
#new_df[is.na(new_df)] <- 0

new_df <- log(new_df[,2:4])
colnames(new_df)[1] <- "Small intestine"
colnames(new_df)[2] <- "Spleen"
colnames(new_df)[3] <- "Duodenum"
transformed <- t(new_df)
d <- dist(as.matrix(transformed))

hc <-hclust(d)
plot(hc)


#removing Freq. from columns names
colnames(new_df) <- gsub("Freq.", "", colnames(new_df))
#removing PTM where both are < 5
#new_df <- new_df[!(new_df$TP_all < 5 & new_df$psm_Crystall < 5),]
#going back to join_df 
summary(new_df)
#frequency_df <- frequency_df[!(frequency_df$Var1==""),]
#making columns as numeric
new_df[2:86] <-lapply(new_df[2:86], as.numeric)

#converting NA to 0
new_df[is.na(new_df)] <- 0
#removing Freq. from columns names
colnames(new_df) <- gsub("Freq.", "", colnames(new_df))

#making row names as sample names
rownames(new_df) <- new_df[,1]
#removing Sample column
new_df$Sample <- NULL

ggplot(data=new_df, aes(x=colnames(new_df), fill=rownames(new_df))) +
  geom_bar(stat = "identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                                                 axis.title = element_text(size=14, face="bold")) +
  labs (x ="PTMs", y = "Frequency (log10)") + ggtitle("title") + scale_y_log10()
#
pClass  <-as.factor(new_df$Sample)


table(pClass)
rownames(new_df) <- new_df[,1]
summary(new_df)




#PCA analysis
Mypca <-prcomp(t(new_df), scale= TRUE)
attributes(Mypca)

s <-summary(Mypca)
str(s)

Xscores   <- Mypca$x

expVar <- s$importance[2,] * 100
expVar

plot(Xscores, xlab="PC1", ylab="PC2", pch=19, cex=1,cex.lab=0.7, cex.axis = 0.7,
     main = "title")

