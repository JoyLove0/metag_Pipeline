######## Functions prepared by Zaid Abdo - 12/2021 - Rota metaomics analysis
### Data preprocessing kraken
data.prepro.ftn = function(file.name){
  d.df = read.table(file = file.name,sep="\t",header = FALSE)
  taxa.names.ls = list()
  row.indx = c()
  rnm.ls = strsplit(d.df[,1],split = "|",fixed=TRUE)
  j = 1
  for(i in 2:length(rnm.ls)){
    #for(i in 2:50){
    if(length(rnm.ls[[i]])>length(rnm.ls[[i-1]])){
      i = i + 1
    }else{
      taxa.names.ls[[j]] = rnm.ls[[i-1]]
      row.indx[j] = i-1
      j = j + 1
    }
  }
  
  # Splitting on __ to identify gaps in the hierarchy of the data
  taxa.split.ls = list()
  for(i in 1:length(taxa.names.ls))
    taxa.split.ls[[i]] = strsplit(taxa.names.ls[[i]],split = "__")
  levels.vt = c("d","k","p","c","o","f","g","s")
  taxa.mt = matrix(nrow=length(taxa.names.ls),ncol=8)
  for(i in 1:length(taxa.names.ls)){
    k = 1
    for(j in 1:length(levels.vt)){
      if(k <= length(taxa.names.ls[[i]])){
        if(levels.vt[j]==taxa.split.ls[[i]][[k]][1]){
          taxa.mt[i,j] = paste(taxa.split.ls[[i]][[k]][1],taxa.split.ls[[i]][[k]][2],sep="__")
          k = k + 1
        }else{
          taxa.mt[i,j] = paste(taxa.split.ls[[i]][[k-1]][1],taxa.split.ls[[i]][[k-1]][2],sep="__")
        }
      }else{
        taxa.mt[i,j] = paste(taxa.split.ls[[i]][[k-1]][1],taxa.split.ls[[i]][[k-1]][2],sep="__")
      }
    }
  }
  taxa.df = data.frame(taxa.mt,d.df[row.indx,2],row.names = d.df[row.indx,1])
  colnames(taxa.df) = c("domain","kingdom","phylum","order","class","family","genus","species","count")
  return(taxa.df)
}
