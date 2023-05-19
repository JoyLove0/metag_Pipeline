# Dependecies 
install.packages("vegan",repos = "http://cran.us.r-project.org")
library(vegan)
source("functions.R")

file.names <- read.table(file = "report.filenames", header = FALSE) #bringing in Kraken2 reports as file.names

# Data cleanup
data.df.ls <- list()
for (i in 1:length(file.names[, 1])) {
  data.df.ls[[i]] <- data.prepro.ftn(file.names[i, ])
}

# Print the length of data.df.ls[[i]] at each iteration
for (i in 1:length(data.df.ls)) {
   print(length(data.df.ls[[i]]))
}

# Data matching
if (length(data.df.ls[[1]]) > length(data.df.ls[[2]])) {
  match.vt <- row.names(data.df.ls[[2]]) %in% row.names(data.df.ls[[1]])
  rnms.vt <- c(row.names(data.df.ls[[1]]), row.names(data.df.ls[[2]])[match.vt == FALSE])
} else {
  match.vt <- row.names(data.df.ls[[1]]) %in% row.names(data.df.ls[[2]])
  rnms.vt <- c(row.names(data.df.ls[[2]]), row.names(data.df.ls[[1]])[match.vt == FALSE])
}

for (i in 3:length(data.df.ls)) {
  if (length(rnms.vt) > length(data.df.ls[[i]])) {
    match.vt <- row.names(data.df.ls[[i]]) %in% rnms.vt
    rnms.vt <- c(rnms.vt, row.names(data.df.ls[[i]])[match.vt == FALSE])
  } else {
    match.vt <- rnms.vt %in% row.names(data.df.ls[[i]])
    rnms.vt <- c(row.names(data.df.ls[[i]]), rnms.vt[match.vt == FALSE])
  }
}

taxa.df <- data.frame(matrix(nrow = length(rnms.vt), ncol = 8))
rownames(taxa.df) <- rnms.vt
colnames(taxa.df) <- c("domain", "kingdom", "phylum", "order", "class", "family", "genus", "species")

data.df <- data.frame(matrix(nrow = length(rnms.vt), ncol = length(data.df.ls)))
rownames(data.df) <- rnms.vt
colnames(data.df) <- sapply(file.names[, 1], function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1])

for (i in 1:length(data.df.ls)) {
  taxa.df[rnms.vt[rnms.vt %in% row.names(data.df.ls[[i]])], ] <- data.df.ls[[i]][rnms.vt[rnms.vt %in% row.names(data.df.ls[[i]])], 1:8]
  data.df[rnms.vt[rnms.vt %in% row.names(data.df.ls[[i]])], i] <- data.df.ls[[i]][rnms.vt[rnms.vt %in% row.names(data.df.ls[[i]])], 9]
}

data.df[is.na(data.df)] <- 0

# Clean-up
rm(data.df.ls, file.names, match.vt, rnms.vt)

# Save the compiled data to a CSV file
write.csv(data.df, "Compiled-Kraken2-Report.csv", row.names = TRUE)

######## Sequence depth
rs.vt = apply(data.df,1,sum) # row sum

pdf(file="depth.pdf",width=4,height=4)
hist(rs.vt,main = "Distribution of Sequencing Depth",xlab = "sequencing depth",ylab = "Number of Samples")
dev.off()
min(rs.vt)
max(rs.vt)

######## Rarifaction curves
pdf(file="rarefaction_curves.pdf", width = 10, height = 10)
rarecurve(x=data.df,col=1:length(rs.vt),xlab = "Number of reads",ylab="OTU")
dev.off()
