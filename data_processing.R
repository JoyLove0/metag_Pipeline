#source("libraries.R")

install.packages("vegan",repos = "http://cran.us.r-project.org")
library(vegan)

source("functions.R")
file.names = read.table(file = "report.filenames",header=FALSE)
# data cleanup
data.df.ls = list()
for(i in 1:length(file.names[,1]))
  data.df.ls[[i]] = data.prepro.ftn(file.names[i,])

# data matching
if(length(data.df.ls[[1]])>length(data.df.ls[[1]])){
  match.vt = row.names(data.df.ls[[2]]) %in% row.names(data.df.ls[[1]])
  rnms.vt = c(row.names(data.df.ls[[1]]),row.names(data.df.ls[[2]])[match.vt==FALSE])
}else{
  match.vt = row.names(data.df.ls[[1]]) %in% row.names(data.df.ls[[2]])
  rnms.vt = c(row.names(data.df.ls[[2]]),row.names(data.df.ls[[1]])[match.vt==FALSE])
}
for(i in 3:length(data.df.ls)){
  if(length(rnms.vt)>length(data.df.ls[[i]])){
    match.vt = row.names(data.df.ls[[i]]) %in% rnms.vt
    rnms.vt = c(rnms.vt,row.names(data.df.ls[[i]])[match.vt==FALSE])
  }else{
    match.vt = rnms.vt %in% row.names(data.df.ls[[i]])
    rnms.vt = c(row.names(data.df.ls[[i]]),rnms.vt[match.vt==FALSE])
  }
}

taxa.df = matrix(nrow = length(rnms.vt),ncol = 8)
taxa.df = data.frame(taxa.df,row.names = rnms.vt)
colnames(taxa.df) = c("domain","kingdom","phylum","order","class","family","genus","species")

data.df = matrix(nrow = length(rnms.vt),ncol = length(data.df.ls))
data.df = data.frame(data.df,row.names = rnms.vt)
df.nms.ls = strsplit(file.names[,1],split = "_",fixed=TRUE)
df.nms = c()
for(i in 1:length(df.nms.ls)) df.nms[i] = df.nms.ls[[i]][1]
colnames(data.df) = df.nms

for(i in 1:length(data.df.ls)){
  taxa.df[rnms.vt[rnms.vt%in%row.names(data.df.ls[[i]])],] = data.df.ls[[i]][rnms.vt[rnms.vt%in%row.names(data.df.ls[[i]])],1:8]
  data.df[rnms.vt[rnms.vt%in%row.names(data.df.ls[[i]])],i] = data.df.ls[[i]][rnms.vt[rnms.vt%in%row.names(data.df.ls[[i]])],9]
}
data.df[is.na(data.df)] = 0
rm(data.df.ls,df.nms.ls,file.names,df.nms,i,match.vt,rnms.vt)

write.csv(data.df,"Compiled-Kraken2-Report.csv")

######## Sequence depth
rs.vt = apply(data.df,1,sum) # row sum (total reads per sample)

pdf(file="depth.pdf",width=4,height=4)
hist(rs.vt,main = "Distribution of Sequencing Depth",xlab = "sequencing depth",ylab = "Number of Samples")
dev.off()
min(rs.vt)
max(rs.vt)

######## Rarifaction curves (we discussed those in class)
#pdf(file="rarefaction_curves.pdf", width = 10, height = 10)
#rarecurve(x=data.df,col=1:length(rs.vt),xlab = "Number of reads",ylab="OTU")
#dev.off()