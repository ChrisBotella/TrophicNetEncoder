repoDir = getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS
source(paste(repoDir,'functions_to_source.R',sep=""))
library(energy)
require(xtable)

masterDir= getwd() #OR REPLACE BY DIRECTORY CONTAINING SCRIPTS
graphsDir = paste(masterDir,'Graphs/',sep="")
ebdDir = paste(masterDir,'Embeddings/',sep="")

# Load property table and embeddings
setwd(ebdDir)
load(file = 'ebds_and_Graphproperties')
# Order of rows in embeddings equal to properties table
orderedNetNames = rownames(ebds[[1]])

#####
# Predictive accuracy (RF OOB accuracy)
# AND variables importance
##### 

caracToPredict = c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")

results = expand.grid(method=c(names(ebds),
                               paste(names(ebds),'RF_rdVars',sep="")),
                      carac=caracToPredict,oob_accuracy=NA)
results$method = as.character(results$method)
results$carac = as.character(results$carac)
results$importanceEntropy = NA
results$Best4SumOverAll = NA

for(i in 1:4){
  eval(parse(text=paste('results$var',i,'=NA',sep="")))
  eval(parse(text=paste('results$scor',i,'=NA',sep="")))
}

for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]]
  
  if(method%in%c('ShortPaths2Vec_lab')){
    ebdSvd = svd(mat)$u[,1:64]
    Xpre = as.matrix(ebdSvd)
  }else{
    Xpre = mat
  }
  
  for(RandomEmbedding in c(F,T)){
    if(!RandomEmbedding){
      X= Xpre
      methodName = as.character(method) 
    }else{
      X = matrix(runif(dim(Xpre)[1]*dim(Xpre)[2],0,1) , dim(Xpre)[1], dim(Xpre)[2])
      methodName = paste(as.character(method),'RF_rdVars',sep="")
    }
    
    for(col in caracToPredict){
      print(caracToPredict[caracToPredict==col])
      Y = factor(exp[,col])
      
      # Select mtry
      nbvars <- seq(1, min(ncol(X) - 1, 13) ,  2) 
      oobsMtry <- sapply(nbvars, function(nbv) {
        RF <- randomForest(Y ~ ., data = X, ntree = 100, mtry = nbv)
        return(RF$err.rate[RF$ntree, "OOB"])})
      
      # RF with optimal mtry
      forest = randomForest(Y ~ ., 
                            data = X,
                            ntree=300,
                            mtry=nbvars[which.min(oobsMtry)],
                            importance=T)
      
      # Store RF OOB accuracy
      results$oob_accuracy[results$method==methodName & results$carac==col] = 1- as.numeric(forest$err.rate[forest$ntree, "OOB"])
      
      if(method=="ShortPaths2Vec_lab"){
        # For SPlab we need the RF computed on original embedding for variables importance
        nbvars <- seq(1, min(ncol(mat) - 1, 13) ,  2) 
        oobsMtry <- sapply(nbvars, function(nbv) {
          RF <- randomForest(Y ~ ., data = mat, ntree = 100, mtry = nbv)
          return(RF$err.rate[RF$ntree, "OOB"])})
        forest = randomForest(Y ~ ., 
                              data = mat,
                              ntree=300,
                              mtry=nbvars[which.min(oobsMtry)],
                              importance=T)}
      
      # Variables importance
      imp_accuracyDecrease = importance(forest,type=1)
      ordered = order(imp_accuracyDecrease,decreasing = T)
      cd = results$method==method & results$carac==col
      results[cd,c(paste('var',1:4,sep=""))] = attr(imp_accuracyDecrease,"dimnames")[[1]][ordered][1:4]
      results[cd,c(paste('scor',1:4,sep=""))] = imp_accuracyDecrease[ordered][1:4]
      
      pis = imp_accuracyDecrease / sum(imp_accuracyDecrease)
      results$importanceEntropy[cd] = - sum( pis * log(pis) )
      results$Best4SumOverAll[cd] = sum(imp_accuracyDecrease[ordered][1:4])/sum(imp_accuracyDecrease[ordered])
    }
  }
  setwd(ebdDir)
  write.table(results,"Segregation_And_VarImp_EbdsAndRandom.csv",sep=";",row.names=F,col.names=T)
}

#####
# Table of predictive accuracy (Table 3)
#####

setwd(ebdDir)
results = read.csv("Segregation_And_VarImp_EbdsAndRandom.csv",sep=";",header=T)

results_tmp = results[regexpr('RF_rdVars',results$method)<0,]

tmp = data.frame(method=unique(results_tmp$method))
tmp$method = as.character(tmp$method)
for(car in unique(results$carac)){
  eval(parse(text=paste('tmp$',car,'=results_tmp$oob_accuracy[results_tmp$cara=="',car,'"]',sep="")))
}
setwd(ebdDir)
write.table(tmp,'segregationTable.csv',sep=";",row.names=F,col.names=T)

rownames(tmp)=tmp$method
tmp = tmp[,colnames(tmp)!="method"]
print(tmp,digits=2)

#####
# Table of R2-Anderson
# on full embeddings
#####

ldists <- lapply(ebds, dist)
props=c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")
res_R2 <- matrix (0, nrow = length(ebds), ncol = length(props))
rownames(res_R2) <- names(ebds)
colnames(res_R2) <- colnames(exp[,props])
for(i in 1:nrow(res_R2)){
  for(j in 1:ncol(res_R2)){
    ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[,props[j]]), permutations = 0)
    res_R2[i, j] <- ado$aov.tab$R2[1]
  }
}

print(res_R2[,colnames(res_R2)!='name'],digits=2)

res_R2= as.data.frame(res_R2)
res_R2$name = rownames(res_R2)

setwd(ebdDir)
write.table(res_R2,'R2anderson_Ebds.csv',sep=";",row.names=F,col.names=T)

#####
# Compute Umap plans 
#####

umaps = list()
for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]]
  
  # reduce dimension to 2 with Umap using default settings
  ebdUmap = as.data.frame(umap(mat,n_neighbors=min(150,dim(mat)[1]))$layout)
  ebdUmap = cbind(ebdUmap,rownames(ebds[[i]]))
  colnames(ebdUmap) = c('umap1','umap2','name')
  umaps[[i]] = ebdUmap
  setwd(ebdDir)
  write.table(ebdUmap,paste("umap_",method,".csv",sep=""),sep=";",row.names=F,col.names=T)
  gc(reset=T)
}
names(umaps) = names(ebds)

umaps = lapply(umaps,function(el)el[,colnames(el)!="name"])

setwd(ebdDir)
save(umaps,file="umaps")

#####
# Table of R2-Anderson
# on 2D plan
#####

setwd(ebdDir)
load(file="umaps")

ldists <- lapply(umaps, dist)
props=c("maxTropLen","TropLens","nchains","omni","generalism","loop","n")
res_R2 <- matrix (0, nrow = length(umaps), ncol = length(props))
rownames(res_R2) <- names(umaps)
colnames(res_R2) <- props
for(i in 1:nrow(res_R2)){
  for(j in 1:ncol(res_R2)){
    ado = vegan::adonis(ldists[[i]] ~ as.factor(exp[,props[j]]), permutations = 0)
    res_R2[i, j] <- ado$aov.tab$R2[1]
  }
}

res_R2= as.data.frame(res_R2)
res_R2$name = rownames(res_R2)

setwd(ebdDir)
write.table(res_R2,'R2anderson_Umaps.csv',sep=";",row.names=F,col.names=T)

#####
# Latex table with all R2 (Table 4)
#####

setwd(ebdDir)
r2_umap = read.csv('R2anderson_Umaps.csv',sep=";",header=T)
r2_ebd = read.csv('R2anderson_Ebds.csv',sep=";",header=T)
rownames(r2_umap) = r2_umap$name
rownames(r2_ebd) = r2_ebd$name

methods = c("Groups2Vec","Metrics","Motifs2Vec","Graph2Vec_dp2","Graph2Vec_lab_dp2",
            "ShortPaths2Vec","ShortPaths2Vec_lab")

tmp = expand.grid(spa=c('embedding','umap'),prop=c('maxTropLen','TropLens','nchains','omni','generalism','loop','n'))
newCols = paste(tmp$spa,'_',tmp$prop,sep="") 

r2_all = matrix(NA,nrow = length(methods),ncol=length(newCols))
for(i in 1:dim(r2_all)[1]){
  for(j in 1:dim(r2_all)[2]){
    if(tmp$spa[j]=="embedding"){
      r2_all[i,j] = r2_ebd[methods[i],tmp$prop[j]]
    }else{
      r2_all[i,j] = r2_umap[methods[i],tmp$prop[j]]
    }
  }
}

r2_all = round(100*r2_all)/100
rownames(r2_all) = methods
colnames(r2_all) = newCols

print(xtable(r2_all))

#####
# Plot Umap plan with colored modalities
# For Figure 2 and 3
##### 

setwd(ebdDir)
load(file="umaps")
oob = read.csv('segregationTable.csv',sep=";",header=T,stringsAsFactors = F)
r2 = read.csv('R2anderson_Umaps.csv',sep=";",header=T,stringsAsFactors = F)

colfunc =colorRampPalette(c("blue", "red"))

df = as.data.frame(expand.grid( method= names(ebds), prop =  c('maxTropLen','nchains','omni','generalism','loop') ))
samp = sample(1:length(orderedNetNames),min(600,dim(ebds[[1]])[1]))
for(i in 1:dim(df)[1]){
  meth = as.character(df$method[i])
  prop =as.character(df$prop[i]) 
  print(meth)
  print(prop)
  
  setwd(ebdDir)
  toPlot = umaps[as.character(df$method[i])][[1]]
  #toPlot = read.csv(paste("umap_",meth,".csv",sep=""),sep=";",header=T,stringsAsFactors = F)
  
  toPlot[,prop] = factor(exp[,prop])
  toPlot$n = factor(exp$n)
  toPlot_tmp=toPlot[samp,]
  
  textToAdd= data.frame(lab=c(paste('oob accur.:',round(100*oob[oob$method==meth,prop])),
                              paste('R2 umap:',round(100*r2[r2$name==meth,prop])/100)),
                        x=rep( min(toPlot_tmp$umap1)+0.12*(max(toPlot_tmp$umap1)-min(toPlot_tmp$umap1)) ,2),
                        y=c( max(toPlot_tmp$umap2)-0.00*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)),max(toPlot_tmp$umap2)-0.03*(max(toPlot_tmp$umap2)-min(toPlot_tmp$umap2)))) 
  
  if(prop=="maxTropLen"){
    plOt = ggplot()+
        geom_point(data=toPlot_tmp,
                   aes(x=toPlot_tmp$umap1,
                       y=toPlot_tmp$umap2,
                       color=toPlot_tmp[,prop],
                       shape=toPlot_tmp$n),size=2)+labs(color=prop)+
        geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=textToAdd$lab))+
        theme_bw()+xlab('umap1')+ylab('umap2')+
        ggtitle(paste('Method:',meth,'; color:',prop))+
        scale_color_manual(values=colfunc(4))+scale_shape_manual(values=c(4,16))
    
  }else{
    plOt = ggplot()+
             geom_point(data=toPlot_tmp,
                        aes(x=toPlot_tmp$umap1,
                            y=toPlot_tmp$umap2,
                            color=toPlot_tmp[,prop],
                            shape=toPlot_tmp$n),size=2)+labs(color=prop)+
             geom_text(data=textToAdd,aes(x=textToAdd$x,y=textToAdd$y,label=textToAdd$lab))+
             theme_bw()+xlab('umap1')+ylab('umap2')+
             ggtitle(paste('Method:',meth,'; color:',prop))+
             scale_color_manual(values=colfunc(2))+scale_shape_manual(values=c(4,16))
  }

  setwd(ebdDir)
  png(paste('Figure_Umap_',meth,'_',prop,'.png',sep=""),
      height=400,
      width=500)
  print(plOt)
  dev.off()
}

#####
# Euclidian distance Scatter plot inter-methods (investigate complementarity)
#####

nSamp = 1000
nSampledDists = 10000

samp = sample(1:dim(ebds[[1]]),nSamp)
NetNames = orderedNetNames[samp]
distances = expand.grid(j=as.integer(1:nSamp),i=as.integer(1:nSamp))
distances = distances[distances$j<distances$i,]
for(i in 1:length(ebds)){
  method = names(ebds)[i]
  print(i)
  print(method)
  mat = ebds[method][[1]][samp,]
  
  # For embeddings, Distance = Euclidean distance between points in the embedding
  Ds = dist(mat,method = "euclidian",diag=F,upper = F)  
  vec = as.numeric(Ds)
  vec = vec/max(vec)
  print('dists computed')
  
  eval(parse(text= paste('distances$',method,'=vec',sep="" ) ))  
}

# We add a row full of 0 to insure that point (0,0) is in all plots windows
toAdd = distances[1,,drop=F]
toAdd[1,] = 0
distances = rbind(toAdd,distances)

sampledDists = c(1,sample(1:dim(distances)[1],nSampledDists))

# Order of methods to plot
refNames = c('Groups2Vec','ShortPaths2Vec_lab','Motifs2Vec','Graph2Vec_lab_dp2','Metrics','ShortPaths2Vec','Graph2Vec_dp2')
toPlotNames = c('Groups2Vec','ShortPaths2Vec_lab','Motifs2Vec','Graph2Vec_lab','Metrics','ShortPaths2Vec','Graph2Vec')

formu = as.formula( paste('~ ', paste(refNames,collapse="+"))) 
print(formu)
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE,breaks="fd")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

# Plot and save multi scatter plot
setwd(ebdDir)
png('methods_distances_pairs_final.png',height=1700,width=1700)
pairs( formu ,data=distances[sampledDists,] , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= toPlotNames,
       pch=3 ,
       col=  alpha(rep('black',dim(distances)[1]),0.05) , 
       cex.labels = 3.5 , 
       cex.axis = 3)
dev.off()

#####
# Distance correlations plot
#####


mcor = matrix(NA,length(refNames),length(refNames))
for(i in 1:length(refNames)){
  for(j in 1:length(refNames)){
    mcor[i,j] = energy::dcor( ebds[refNames[i]][[1]] , ebds[refNames[j]][[1]])
  }
}
rownames(mcor) = toPlotNames
colnames(mcor) = toPlotNames

setwd(ebdDir)
saveRDS(mcor,'mat_dCor')

#print(mcor,digits=2)

setwd(ebdDir)
png('methods_distances_correlations_final.png',height=2000,width=2000)
print(corrplot(mcor, type="upper", tl.col="white", tl.srt=45, diag=F,cl.cex=3,cl.lim=c(0.25,1),addCoef.col = "grey80",number.cex=4))
dev.off()

# Old agreement metric
if(F){
  library(corrplot)
  
  mcor = matrix(NA,length(toPlot),length(toPlot))
  for(i in 1:length(toPlot)){
    for(j in 1:length(toPlot)){
      mcor[i,j] = cor(distances[,as.character(names(ebds)[toPlot[i]])],distances[,as.character(names(ebds)[toPlot[j]])])
    }
  }
  rownames(mcor) = names(ebds)[toPlot]
  colnames(mcor) = names(ebds)[toPlot]
  
  print(mcor,digits=2)
  
  setwd(ebdDir)
  png('methods_distances_pairs_correlations_OLD.png',height=2000,width=2000)
  print(corrplot(mcor[reord,reord], type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=3))
  dev.off()
  
}
