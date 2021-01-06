#####
# Table of R2-loss in embedding
#####

Ns = c(60,120)
resList = list()
for(N in Ns){
  SAMPLE = exp$n==N
  ldists <- lapply(1:length(ebds),function(i) dist(ebds[[i]][SAMPLE,])  )
  props=c("maxTropLen","TropLens","nchains","omni","generalism","loop")
  res_R2 <- matrix (0, nrow = length(ebds), ncol = length(props))
  rownames(res_R2) <- names(ebds)
  colnames(res_R2) <- colnames(exp[,props])
  for(i in 1:nrow(res_R2)){
    for(j in 1:ncol(res_R2)){
      ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[SAMPLE,props[j]]), permutations = 0)
      res_R2[i, j] <- ado$aov.tab$R2[1]
    }
  }
  resList[[which(Ns==N)]] = res_R2
}

results = (resList[[1]] + resList[[2]])/2
results= as.data.frame(results)
results$name = rownames(results)

setwd(ebdDir)
write.table(results,'MeanR2anderson_Ebds.csv',sep=";",row.names=F,col.names=T)


# Print R2-loss table
setwd(ebdDir)
R2mean = read.csv('MeanR2anderson_Ebds.csv',sep=";",header=T)
R2mixed = read.csv('R2anderson_Ebds.csv',sep=";",header=T) # Output from FullAnalysisPipeline.R

meths = as.character(R2mixed$name)

R2mean= R2mean[,colnames(R2mean)!="name"]
R2mixed = R2mixed[,!colnames(R2mixed)%in%c('name','n')]

print('R2_loss full embedding')
res = (R2mean - R2mixed)/R2mean
rownames(res)=meths
print( round(100*res),digits=2)

#####
# Table of R2-loss on Umap
#####

Ns = c(60,120)
resList = list()
for(N in Ns){
  SAMPLE = exp$n==N
  ldists <- lapply(1:length(umaps),function(i) dist(umaps[[i]][SAMPLE,])  )
  props=c("maxTropLen","TropLens","nchains","omni","generalism","loop")
  res_R2 <- matrix (0, nrow = length(umaps), ncol = length(props))
  rownames(res_R2) <- names(umaps)
  colnames(res_R2) <- colnames(exp[,props])
  for(i in 1:nrow(res_R2)){
    for(j in 1:ncol(res_R2)){
      ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[SAMPLE,props[j]]), permutations = 0)
      res_R2[i, j] <- ado$aov.tab$R2[1]
    }
  }
  resList[[which(Ns==N)]] = res_R2
}

results = (resList[[1]] + resList[[2]])/2
results= as.data.frame(results)
results$name = rownames(results)

setwd(ebdDir)
write.table(results,'MeanR2anderson_Umaps.csv',sep=";",row.names=F,col.names=T)

# Print R2-loss in Umap plan
setwd(ebdDir)
R2mean = read.csv('MeanR2anderson_slicedEbds.csv',sep=";",header=T)
R2mixed = read.csv('R2anderson_slicedEbds.csv',sep=";",header=T)
meths = as.character(R2mixed$name)
R2mean= R2mean[,colnames(R2mean)!="name"]
R2mixed = R2mixed[,!colnames(R2mixed)%in%c('name','n')]
res = (R2mean - R2mixed)/R2mean
rownames(res)=meths
print( round(100*res),digits=2)

