install.packages("ggplot2")
library(ggplot2)

df <- read.csv("elapsed_time.csv", sep = ";")


binarysearch <- subset(df, df$Step < 3) #obtaining different rows where "Step" is less than 3

ggplot(binarysearch, aes(x=Step, y=Cumulative_time, group=Algorithm, color=Algorithm))+
geom_line(size=0.3) + geom_point(size= 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9),
labels = c("Start","reading file", "Binary search (1 PTM)", "loop search control (2PTM)",
"loop search control (3PTM)", "join 2 and 3 PTM search", "loop search control (4PTM)",
 "join 4 and 5 natural PTM search", "removing search control", "join 4 and 5 PTM search")) +
labs(y= "Time (seconds)", x = "Improvement steps") + ggtitle("Execution time") + theme(plot.title = element_text(hjust = 0.5))



#difference of time from moving in steps
ggplot(binarysearch, aes(x=Step, y=Dif, group=Algorithm, color=Algorithm))+
geom_line(size=0.3) + geom_point(size= 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
 scale_x_continuous(breaks = c(1,2),
labels = c("reading file", "Binary search (1 PTM)")) +
labs(y= "Time (seconds)", x = "Improvement steps") + ggtitle("Execution time") + theme(plot.title = element_text(hjust = 0.5))
