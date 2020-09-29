install.packages("ggplot2")
library(ggplot2)

df <- read.csv("/Users/pedrocardoso/Documents/PTMs_project/elapsed_time.csv", sep = ";")


binarysearch <- subset(df, df$Step > 0 & df$Step <11) #obtaining different rows where "Step" is less than 3

ggplot(binarysearch, aes(x=Step, y=Cumulative_time, group=Algorithm, color=Algorithm))+
geom_line(size=0.3) + geom_point(size= 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 scale_x_continuous(breaks = c(1,2,3,4,5,6),
labels = c("reading file", "Binary search (1 PTM)", "loop search control (2PTM)",
"loop search control (3PTM)", "join 2 and 3 PTM search", "loop search control (4PTM)", )) +
labs(y= "Time (seconds)", x = "Improvement steps") + ggtitle("Execution time") +
theme(axis.title = element_text(size=14, face="bold"),
      axis.text = element_text(size=12),
      legend.text = element_text(size=12),
      legend.title = element_text(size=14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5))



#difference of time from moving in steps
ggplot(binarysearch, aes(x=Step, y=Dif, group=Algorithm, color=Algorithm))+
geom_line(size=0.3) + geom_point(size= 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 scale_x_continuous(breaks = c(1,2),
labels = c("reading file", "Binary search (1 PTM)")) +
labs(y= "Time (seconds)", x = "Improvement steps") + ggtitle("Execution time") + theme(plot.title = element_text(hjust = 0.5))
