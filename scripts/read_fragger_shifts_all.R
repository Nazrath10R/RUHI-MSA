
L1_digested <- read.table("~/bessantlab/naz/biological_interpretation/sequences/LINE_1_digested.tsv",
                          header = TRUE)

path <- "~/Desktop/ptm/ptm/ptm/"

setwd("~/bessantlab/naz/fragpipe/results/deamidation/")

df <- data.frame("Sample"=NA, "Spectrum"=NA, "Peptide"=NA, "Modified.Peptide"=NA, "Charge"=NA, 
                 "Retention"=NA, "Calculated.M.Z"=NA, "Observed.M.Z"=NA, "Original.Delta.Mass"=NA, 
                 "Adjusted.Delta.Mass"=NA, "Experimental.Mass"=NA, "Peptide.Mass"=NA, "Expectation"=NA, 
                 "Hyperscore"=NA, "Nextscore"=NA, "PeptideProphet.Probability"=NA, "Intensity"=NA, 
                 "Assigned.Modifications"=NA, "Observed.Modifications"=NA, "Number.of.Phospho.Sites"=NA, 
                 "Phospho.Site.Localization"=NA, "Is.Unique"=NA, "Protein"=NA, "Protein.ID"=NA, "Entry.Name"=NA, 
                 "Gene"=NA, "Protein.Description"=NA, "Mapped.Proteins"=NA)



dirs <- list.dirs(recursive = FALSE)
dirs <- dirs[-1]


for(x in 1:length(dirs)) {
  
  psm_table <- read.table(paste(dirs[x], "/psm.tsv", sep = ""), 
                          sep = "\t", header = TRUE, quote="")
  
  # L1_peptides <- intersect(psm_table$Peptide, L1_digested$Peptide)
  # L1_table_subset <- subset(psm_table, Peptide %in% L1_peptides)
  # L1_table_subset <- cbind("Sample"=rep(sample, nrow(L1_table_subset)),L1_table_subset)
  
  sample <- unlist(strsplit(dirs[x], split = "./"))[2]
  
  table_subset <- cbind("Sample"=rep(sample, nrow(psm_table)),psm_table)
  
  df <- rbind(df, table_subset)
  
}

df <- df[-1,]


df <- df[which(df$PeptideProphet.Probability > 0.99),]
df_mass_shifted <- df[which(df$Assigned.Modifications==""),]



names <- df_mass_shifted$Peptide

numbers <- format(as.numeric(as.character(df_mass_shifted$Adjusted.Delta.Mass)), nsmall = 4, scientific = FALSE)
numbers <- gsub(" ", "",numbers)

masses <- format(df_mass_shifted$Observed.M.Z, nsmall = 6, scientific = FALSE)
masses <- gsub(" ", "",masses)

output_table <- c(rbind(names, numbers, masses))



sample_name <- "all_fragger_pp"

dir.create(paste(path, "results/", sample_name, sep = ""))

output_name <- paste(path, "input/", sample_name, ".txt", sep = "")

# df_mass_shifted

write.table(output_table, output_name, quote = FALSE, row.names = FALSE,  col.names = FALSE)







# for(x in 1:length(unique(df$Sample))) {
#   sample_x <- df[df$Sample==unique(df$Sample)[x],]
#   print(unique(sample_x$Adjusted.Delta.Mass - sample_x$Original.Delta.Mass))
# 
# }

































